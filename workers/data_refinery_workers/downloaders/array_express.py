from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import zipfile
import time
from typing import List
from contextlib import closing

from data_refinery_common.models import (
    DownloaderJob,
    Experiment,
    Sample,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    DownloaderJobOriginalFileAssociation
)
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import (
    get_env_variable,
    get_readable_affymetrix_names,
    get_supported_microarray_platforms
)
from data_refinery_common import microarray


logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    """ Download a file from ArrayExpress via FTP. There is no Aspera endpoint
    which I can find. """
    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path,
                     downloader_job=job.id)
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url, timeout=60)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logger.exception("Exception caught while downloading file.",
                         downloader_job=job.id)
        job.failure_reason = "Exception caught while downloading file"
        raise
    finally:
        target_file.close()


def _extract_files(file_path: str, accession_code: str, job: DownloaderJob) -> List[str]:
    """Extract zip and return a list of the raw files.
    """

    logger.debug("Extracting %s!", file_path, file_path=file_path, downloader_job=job.id)

    try:
        # This is technically an unsafe operation.
        # However, we're trusting AE as a data source.
        zip_ref = zipfile.ZipFile(file_path, "r")
        abs_with_code_raw = LOCAL_ROOT_DIR + '/' + accession_code + '/raw/'
        zip_ref.extractall(abs_with_code_raw)
        # Other zips for this same accession will go into this
        # directory too, so look at what's in the zip file rather than
        # what's in the directory it's being extracted to.
        files_in_zip = zip_ref.namelist()
        zip_ref.close()

        # os.abspath doesn't do what I thought it does, hency this monstrocity.
        files = [{'absolute_path': abs_with_code_raw + f, 'filename': f} for f in files_in_zip]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), str(file_path)
        logger.exception(reason, downloader_job=job.id)
        job.failure_reason = reason
        raise

    return files


def download_array_express(job_id: int) -> None:
    """The main function for the Array Express Downloader.

    Downloads a single zip file containing the .PCL files representing
    samples relating to a single experiement stored in
    ArrayExpress.
    """
    job = utils.start_job(job_id)
    success = True

    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
    # AE will have multiple files per DownloaderJob, but they are all
    # pieces of the same zip file so they're all referencing the same
    # URL.
    original_file = file_assocs[0].original_file
    url = original_file.source_url
    accession_code = job.accession_code

    # First, get all the unique sample archive URLs.
    # There may be more than one!
    # Then, unpack all the ones downloaded.
    # Then create processor jobs!

    og_files = []
    # The files for all of the samples are
    # contained within the same zip file. Therefore only
    # download the one.
    os.makedirs(LOCAL_ROOT_DIR + '/' + accession_code, exist_ok=True)

    # Add a timestamp in milliseconds to filename to prevent multiple jobs from using the same file.
    filename = url.split('/')[-1] + "." + str(int(time.time() * 1000))
    dl_file_path = LOCAL_ROOT_DIR + '/' + accession_code + '/' + filename + ".zip"
    _download_file(url, dl_file_path, job)

    extracted_files = _extract_files(dl_file_path, accession_code, job)

    for og_file in extracted_files:
        # TODO: We _should_ be able to use GET here - anything more than 1 sample per
        # filename is a problem. However, I need to know more about the naming convention.
        try:
            original_file = OriginalFile.objects.filter(
                source_filename=og_file['filename']).order_by('created_at')[0]
            original_file.is_downloaded = True
            original_file.is_archive = False
            original_file.absolute_file_path = og_file['absolute_path']
            original_file.filename = og_file['absolute_path'].split('/')[-1]
            original_file.calculate_size()
            original_file.save()
            original_file.calculate_sha1()
            og_files.append(original_file)
        except Exception:
            # TODO - is this worth failing a job for?
            logger.warn("Found a file we didn't have an OriginalFile for! Why did this happen?: "
                        + og_file['filename'],
                        downloader_job=job_id)
            continue

        sample_objects = Sample.objects.filter(originalfile=original_file).order_by('created_at')
        if sample_objects.count() > 1:
            logger.warn("Found an Array Express OriginalFile with more than one sample: %s",
                        filename,
                        downloader_job=job_id)

        # If the file is a .CEL file, it is the ultimate
        # source of truth about the sample's platform.
        sample_object = sample_objects[0]
        if og_file["filename"].upper()[-4:] == ".CEL" and sample_object.has_raw:
            platform_accession_code = "UNSUPPORTED"
            try:
                cel_file_platform = microarray.get_platform_from_CEL(
                    original_file.absolute_file_path)

                for platform in get_supported_microarray_platforms():
                    if platform["platform_accession"] == cel_file_platform:
                        platform_accession_code = platform["platform_accession"]
            except Exception:
                platform_accession_code = "UNDETERMINABLE"
                logger.warn("Unable to determine platform from CEL file: "
                            + original_file.absolute_file_path,
                            downloader_job=job_id)

            if platform_accession_code == "UNSUPPORTED":
                logger.error("Found a raw .CEL file with an unsupported platform!",
                             file_name=original_file.absolute_file_path,
                             sample=sample_object.id,
                             downloader_job=job_id)
                job.failure_reason = ("Found a raw .CEL file with an unsupported platform: "
                                      + original_file.absolute_file_path)
                success = False
            elif platform_accession_code == "UNDETERMINABLE":
                # If we cannot determine the platform from the
                # .CEL file, the platform discovered via metadata
                # may be correct so just leave it be.
                pass
            else:
                # We determined the file was collected with a supported Affymetrix platform.
                sample_object.platform_accession_code = platform_accession_code
                sample_object.platform_name = get_readable_affymetrix_names()[
                    platform_accession_code]

            # However, if the filename contains '.CEL' we know
            # it's an Affymetrix Microarray
            sample_object.technology = "MICROARRAY"
            sample_object.manufacterer = "AFFYMETRIX"
            sample_object.save()

    if success:
        logger.debug("File downloaded and extracted successfully.",
                     url,
                     downloader_job=job_id)

        utils.create_processor_jobs_for_original_files(og_files, job)

    utils.end_downloader_job(job, success)

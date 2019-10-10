import os
import shutil
import time
import urllib.request
import zipfile

from contextlib import closing
from typing import List

from data_refinery_common import microarray
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    Sample,
)
from data_refinery_common.utils import (
    get_env_variable,
    get_readable_affymetrix_names,
    get_supported_microarray_platforms,
)
from data_refinery_workers.downloaders import utils
from data_refinery_common.job_management import create_processor_jobs_for_original_files


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
    abs_with_code_raw = LOCAL_ROOT_DIR + '/' + accession_code + '/raw/'

    try:
        # This is technically an unsafe operation.
        # However, we're trusting AE as a data source.
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            zip_ref.extractall(abs_with_code_raw)
            # Other zips for this same accession will go into this
            # directory too, so look at what's in the zip file rather than
            # what's in the directory it's being extracted to.
            files_in_zip = zip_ref.namelist()

        return [{'absolute_path': abs_with_code_raw + f, 'filename': f} for f in files_in_zip]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), str(file_path)
        logger.exception(reason, downloader_job=job.id)
        job.failure_reason = reason
        raise


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

    unprocessed_original_files = []
    # The files for all of the samples are
    # contained within the same zip file. Therefore only
    # download the one.
    os.makedirs(LOCAL_ROOT_DIR + '/' + accession_code, exist_ok=True)

    # Add a timestamp in milliseconds to filename to prevent multiple jobs from using the same file.
    filename = url.split('/')[-1] + "." + str(int(time.time() * 1000))
    dl_file_path = LOCAL_ROOT_DIR + '/' + accession_code + '/' + filename + ".zip"
    _download_file(url, dl_file_path, job)

    extracted_files = _extract_files(dl_file_path, accession_code, job)
    os.remove(dl_file_path) # remove zip file

    for extracted_file in extracted_files:
        try:
            original_file = OriginalFile.objects.get(
                source_filename=extracted_file['filename'], source_url=original_file.source_url)
            # Sometimes a file needs to be redownloaded and processed,
            # but if a file is part of an archive and then we requeue
            # all of the files to be processed, we might end up
            # reprocessing files. Therefore make sure that they
            # haven't actually been processed before marking them as
            # downloaded and queuing processor jobs.
            if original_file.needs_processing():
                original_file.set_downloaded(extracted_file['absolute_path'])
                unprocessed_original_files.append(original_file)
        except Exception:
            # The suspicion is that there are extra files related to
            # another experiment, that we don't want associated with
            # this one.
            logger.debug("Found a file we didn't have an OriginalFile for! Why did this happen?",
                        file_name=original_file.filename,
                        downloader_job=job_id)
            os.remove(extracted_file["absolute_path"])
            continue

        sample_objects = original_file.samples.order_by('created_at')
        if sample_objects.count() > 1:
            logger.warn("Found an Array Express OriginalFile with more than one sample",
                        original_file = original_file,
                        downloader_job=job_id)

        # If the file is a .CEL file, it is the ultimate
        # source of truth about the sample's platform.
        sample_object = sample_objects.first()
        if extracted_file["filename"].upper()[-4:] == ".CEL" and sample_object.has_raw:
            cel_file_platform = None
            platform_accession_code = "UNSUPPORTED"
            try:
                cel_file_platform = microarray.get_platform_from_CEL(
                    original_file.absolute_file_path)

                for platform in get_supported_microarray_platforms():
                    if platform["platform_accession"] == cel_file_platform:
                        platform_accession_code = platform["platform_accession"]
            except Exception as e:
                platform_accession_code = "UNDETERMINABLE"
                logger.warn("Unable to determine platform from CEL file: "
                            + str(original_file.absolute_file_path),
                            downloader_job=job_id)
            if platform_accession_code == "UNSUPPORTED":
                logger.error("Found a raw .CEL file with an unsupported platform!",
                             file_name=original_file.absolute_file_path,
                             sample=sample_object.id,
                             downloader_job=job_id,
                             cel_file_platform=cel_file_platform)
                job.failure_reason = ("Found a raw .CEL file with an unsupported platform: "
                                      + original_file.absolute_file_path + " ("
                                      + str(cel_file_platform) + ")")
                job.no_retry = True
                success = False

                # The file is unsupported, delete it!
                original_file.delete_local_file()
                original_file.delete()
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
                     url=url,
                     downloader_job=job_id)

        create_processor_jobs_for_original_files(unprocessed_original_files, job)

    utils.end_downloader_job(job, success)

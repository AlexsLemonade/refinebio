from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import gzip
import tarfile
from typing import List
from contextlib import closing

from data_refinery_common.models import (
    DownloaderJob,
    Experiment, 
    Sample, 
    ExperimentAnnotation, 
    ExperimentSampleAssociation, 
    OriginalFile, 
    DownloaderJobOriginalFileAssociation,
    OriginalFileSampleAssociation
)
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable


logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256

def get_miniml_url(self, experiment_accession_code):
    """ """
    geo = experiment_accession_code.upper()
    geotype = geo[:3]
    range_subdir = sub(r"\d{1,3}$", "nnn", geo)

    miniml_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
              "{root}/{range_subdir}/{record}/miniml/{record_file}")
    miniml_url = miniml_url_template.format(root="series",
                        range_subdir=range_subdir,
                        record=geo,
                        record_file="%s_family.xml.tgz" % geo)

    return raw_url

def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    """ Download a file from ArrayExpress via FTP. There is no Aspera endpoint
    which I can find. """
    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path,
                     downloader_job=job.id)
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logger.exception("Exception caught while downloading file.",
                         downloader_job=job.id)
        job.failure_reason = "Exception caught while downloading file"
        raise
    finally:
        target_file.close()

def _extract_tar(file_path: str, accession_code: str) -> List[str]:
    """Extract tar and return a list of the raw files.
    """

    logger.debug("Extracting %s!", file_path, file_path=file_path)

    try:
        # This is technically an unsafe operation.
        # However, we're trusting GEO as a data source.
        zip_ref = tarfile.TarFile(file_path, "r")
        abs_with_code_raw = LOCAL_ROOT_DIR + '/' + accession_code + '/raw/'
        zip_ref.extractall(abs_with_code_raw)
        zip_ref.close()

        # os.abspath doesn't do what I thought it does, hency this monstrocity.
        files = [{'absolute_path': abs_with_code_raw + f, 'filename': f} for f in os.listdir(abs_with_code_raw)]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), file_path
        logger.exception(reason)
        raise

    return files

def _extract_tgz(file_path: str, accession_code: str) -> List[str]:
    """Extract tgz and return a list of the raw files.
    """

    logger.debug("Extracting %s!", file_path, file_path=file_path)

    try:
        extracted_filepath = file_path.replace('.tgz', '.tar')
        with gzip.open(file_path, 'rb') as f_in:
            with open(extracted_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        zip_ref = tarfile.TarFile(extracted_filepath, "r")
        abs_with_code_raw = LOCAL_ROOT_DIR + '/' + accession_code + '/raw/'
        zip_ref.extractall(abs_with_code_raw)
        zip_ref.close()

        files = [{'absolute_path': abs_with_code_raw + f, 'filename': f} for f in os.listdir(abs_with_code_raw)]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), file_path
        logger.exception(reason)
        raise

    return files

def _extract_gz(file_path: str, accession_code: str) -> List[str]:
    """Extract gz and return a list of the raw files.
    """

    logger.debug("Extracting %s!", file_path, file_path=file_path)

    try:

        extracted_filepath = file_path.replace('.gz', '')
        with gzip.open(file_path, 'rb') as f_in:
            with open(extracted_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        # os.abspath doesn't do what I thought it does, hency this monstrocity.
        files = [{'absolute_path': extracted_filepath, 'filename': extracted_filepath.split('/')[-1]}]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), file_path
        logger.exception(reason)
        raise

    return files

def download_geo(job_id: int) -> None:
    """The main function for the GEO Downloader.

    Downloads a single tar file containing the files representing
    samples relating to a single experiement stored in
    GEO.
    """
    job = utils.start_job(job_id)

    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
    original_file = file_assocs[0].original_file # GEO should never have more than one zip, but we can iterate here if we discover this is false.
    url = original_file.source_url
    accession_code = job.accession_code

    # First, get all the unique sample archive URLs.
    # There may be more than one!
    # Then, unpack all the ones downloaded.
    # Then create processor jobs!


    # The files for all of the samples are
    # contained within the same zip file. Therefore only
    # download the one.
    os.makedirs(LOCAL_ROOT_DIR + '/' + accession_code, exist_ok=True)
    dl_file_path = LOCAL_ROOT_DIR + '/' + accession_code + '/' + url.split('/')[-1]
    _download_file(url, dl_file_path, job)
    original_file.is_downloaded = True
    original_file.save()

    # This files are tarred, and also subsequently gzipped
    if '.tar' in dl_file_path:
        extracted_files = _extract_tar(dl_file_path, accession_code)

        unpacked_sample_files = []
        for og_file in extracted_files:

            filename = og_file['filename']
            sample_id = filename.split('_')[0]

            try:
                # Files from the GEO supplemental file are gzipped inside of the tarball. Great!
                archive_file = OriginalFile.objects.get(source_filename__contains=sample_id)
                archive_file.is_downloaded=True
                archive_file.is_archive=True
                archive_file.absolute_file_path = og_file['absolute_path']
                archive_file.calculate_size()
                archive_file.calculate_sha1()
                archive_file.save()

                extracted_subfile = _extract_gz(og_file['absolute_path'], accession_code)
                actual_file = OriginalFile()
                actual_file.is_downloaded=True
                actual_file.is_archive=False
                actual_file.absolute_file_path = extracted_subfile[0]['absolute_path']
                actual_file.filename = extracted_subfile[0]['filename']
                actual_file.calculate_size()
                actual_file.calculate_sha1()
                actual_file.has_raw = True
                actual_file.sample = archive_file.sample
                actual_file.save()

                # Question - do we want to delete this extracted archive file?
                # archive_file.delete()

                unpacked_sample_files.append(actual_file)
            except Exception as e:
                # TODO - is this worth failing a job for?
                logger.warn("Found a file we didn't have an OriginalFile for! Why did this happen?: " + og_file['filename'])

    # This is a .tgz file.
    elif '.tgz' in dl_file_path:
        extracted_files = _extract_tgz(dl_file_path, accession_code)
        unpacked_sample_files = []
        for og_file in extracted_files:

            if '.txt' in og_file['filename']:
                try:
                    gsm_id = og_file['filename'].split('-')[0]
                    sample = Sample.objects.get(accession_code=gsm_id)
                except Exception as e:
                    print(e)
                    print(og_file)
                    print
                    continue

                actual_file = OriginalFile()
                actual_file.is_downloaded=True
                actual_file.is_archive=False
                actual_file.absolute_file_path = og_file['absolute_path']
                actual_file.filename = og_file['filename']
                actual_file.calculate_size()
                actual_file.calculate_sha1()
                actual_file.has_raw = True
                actual_file.save()

                original_file_sample_association = OriginalFileSampleAssociation()
                original_file_sample_association.sample = sample
                original_file_sample_association.original_file = actual_file
                original_file_sample_association.save()

                unpacked_sample_files.append(actual_file)

            else:
                print(og_file['filename'])

    # These files are only gzipped.
    else:
        extracted_files = _extract_gz(dl_file_path, accession_code)
        unpacked_sample_files = []
        for og_file in extracted_files:

            filename = og_file['filename']
            sample_id = filename.split('.')[0]

            try:
                # The archive we downloaded
                archive_file = OriginalFile.objects.get(source_filename__contains=filename)
                archive_file.is_downloaded=True
                archive_file.is_archive=True
                archive_file.absolute_file_path = dl_file_path
                archive_file.calculate_size()
                archive_file.calculate_sha1()
                archive_file.save()

                actual_file = OriginalFile()
                actual_file.is_downloaded=True
                actual_file.is_archive=False
                actual_file.absolute_file_path = og_file['absolute_path']
                actual_file.filename = og_file['filename']
                actual_file.calculate_size()
                actual_file.calculate_sha1()
                actual_file.has_raw = True
                actual_file.sample = archive_file.sample
                actual_file.save()

                # Question - do we want to delete this extracted archive file?
                # archive_file.delete()

                unpacked_sample_files.append(actual_file)
            except Exception as e:
                # TODO - is this worth failing a job for?
                logger.warn("Found a file we didn't have an OriginalFile for! Why did this happen?: " + og_file['filename'])

    if len(unpacked_sample_files) > 0:
        success = True
        logger.debug("File downloaded and extracted successfully.",
                     url,
                     downloader_job=job_id)
    else:
        success = False
        logger.debug("Unable to extract any files.",
                     url,
                     downloader_job=job_id)

    utils.end_downloader_job(job, success)

    if success:

        # We're trying to detect technology type here.
        # It may make more sense to try to make this into a higher level Sample property.
        annotations = actual_file.samples.first().sampleannotation_set.all()[0]

        # XXX: Make sure this still works if we get arrays back. Should check for the presence of the 'Cy5' string in label_ch2.
        if ('Agilent' in annotations.data.get('label_protocol_ch1', "")) and ('Agilent' in annotations.data.get('label_protocol_ch2', "")):
            utils.create_processor_jobs_for_original_files(unpacked_sample_files, pipeline="AGILENT_TWOCOLOR_TO_PCL")
        if ('Illumina' in annotations.data.get('label_protocol_ch1', "")):
            utils.create_processor_jobs_for_original_files(unpacked_sample_files, pipeline="ILLUMINA_TO_PCL")
        else:
            utils.create_processor_jobs_for_original_files(unpacked_sample_files)

import gzip
import os
import shutil
import subprocess
import tarfile
import time
import urllib.request

from contextlib import closing
from typing import List, Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    Sample,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.downloaders import utils


logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _download_file(download_url: str, file_path: str, job: DownloaderJob, force_ftp=False) -> None:
    """ Download a file from GEO via Aspera unless `force_ftp` is True
    """

    # Ensure directory exists
    os.makedirs(file_path.rsplit('/', 1)[0], exist_ok=True)

    if not force_ftp:
        return _download_file_aspera(download_url=download_url, downloader_job=job, target_file_path=file_path)
    else:
        try:
            logger.debug("Downloading file from %s to %s.",
                         download_url,
                         file_path,
                         downloader_job=job.id)

            # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
            urllib.request.urlcleanup()

            target_file = open(file_path, "wb")
            with closing(urllib.request.urlopen(download_url)) as request:
                shutil.copyfileobj(request, target_file, CHUNK_SIZE)

            urllib.request.urlcleanup()
        except Exception:
            logger.exception("Exception caught while downloading file.",
                             downloader_job=job.id)
            job.failure_reason = "Exception caught while downloading file"
            raise
        finally:
            target_file.close()

        return True


def _download_file_aspera(download_url: str,
                          downloader_job: DownloaderJob,
                          target_file_path: str,
                          attempt=0) -> bool:
    """ Download a file to a location using Aspera by shelling out to the `ascp` client. """

    try:
        logger.debug("Downloading file from %s to %s via Aspera.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        ascp = ".aspera/cli/bin/ascp"
        key = ".aspera/cli/etc/asperaweb_id_dsa.openssh"
        url = download_url
        user = "anonftp"
        ftp = "ftp-trace.ncbi.nlm.nih.gov"
        if url.startswith("ftp://"):
            url = url.replace("ftp://", "")
        url = url.replace(ftp, "").replace('ftp.ncbi.nlm.nih.gov', '')

        # Resume level 1, use encryption, unlimited speed
        command_str = "{} -i {} -k1 -T {}@{}:{} {}".format(ascp, key, user, ftp, url, target_file_path)
        formatted_command = command_str.format(src=download_url,
                                               dest=target_file_path)
        completed_command = subprocess.run(formatted_command.split(),
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

        # Something went wrong! Else, just fall through to returning True.
        if completed_command.returncode != 0:

            stderr = completed_command.stderr.decode().strip()
            logger.debug("Shell call of `%s` to ascp failed with error message: %s",
                         formatted_command,
                         stderr,
                         downloader_job=downloader_job.id)

            # Sometimes, GEO fails mysteriously.
            # Wait a few minutes and try again.
            if attempt >= 5:
                downloader_job.failure_reason = stderr
                logger.error("All attempts to download accession via ascp failed: %s\nCommand was: %s",
                             stderr,
                             formatted_command,
                             downloader_job=downloader_job.id)
                return False
            else:
                time.sleep(30)
                return _download_file_aspera(download_url,
                                             downloader_job,
                                             target_file_path,
                                             attempt + 1
                                             )
    except Exception:
        logger.exception("Exception caught while downloading file from the URL via Aspera: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "file from the URL via Aspera: {}").format(download_url)
        return False

    # If Aspera has given a zero-byte file for some reason, let's back off and retry.
    if os.path.getsize(target_file_path) < 1:
        os.remove(target_file_path)
        if attempt > 5:
            downloader_job.failure_reason = "Got zero byte file from aspera after 5 attempts."
            return False

        logger.error("Got zero byte ascp download for target, retrying.",
                     target_url=download_url,
                     downloader_job=downloader_job.id)
        time.sleep(10)
        return _download_file_aspera(download_url,
                                     downloader_job,
                                     target_file_path,
                                     attempt + 1
                                     )
    return True


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
        files = [{'absolute_path': abs_with_code_raw + f, 'filename': f}
                 for f in os.listdir(abs_with_code_raw)]

    except Exception as e:
        logger.exception("While extracting %s caught exception %s",
                         file_path,
                         str(e),
                         accession_code=accession_code,
                         file_path=file_path)
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

        files = [{'absolute_path': abs_with_code_raw + f, 'filename': f}
                 for f in os.listdir(abs_with_code_raw)]

    except Exception as e:
        reason = "Exception %s caught while extracting %s", str(e), file_path
        logger.exception(reason, accession_code=accession_code, file_path=file_path)
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

        files = [{'absolute_path': extracted_filepath,
                  'filename': extracted_filepath.rsplit('/', 1)[1]
                  }]

    except Exception as e:
        logger.exception("While extracting %s caught exception %s",
                         file_path,
                         str(e),
                         accession_code=accession_code,
                         file_path=file_path)
        raise

    return files


def _get_actual_file_if_queueable(
        extracted_subfile: Dict,
        original_file: OriginalFile,
        samples: List[Sample]) -> OriginalFile:
    """Returns the actual file from the archive if it should be queued.

    If the file has been processed or has an unstarted DownloaderJob,
    None will be returned.

    `extracted_subfile` should be a Dict containing metadata about the
    file that was extracted from an archive.

    `original_file` should be the file associated with the CURRENT
    DownloaderJob.

    `samples` are the samples that the actual file should be associated
    with if it has to be created.
    """
    # Check to see if we've made this original file before:
    potential_existing_files = OriginalFile.objects.filter(
        source_filename=original_file.source_filename,
        filename=extracted_subfile['filename'],
        is_archive=False
    )
    if potential_existing_files.count() > 0:
        # We've already created this record, let's see if we actually
        # needed to download it or if we just got it because we needed
        # a file in the same archive.
        actual_file = potential_existing_files[0]
        if actual_file.needs_downloading():
            return actual_file
        else:
            return None
    else:
        actual_file = OriginalFile()
        actual_file.is_downloaded = True
        actual_file.is_archive = False
        actual_file.absolute_file_path = extracted_subfile['absolute_path']
        actual_file.filename = extracted_subfile['filename']
        actual_file.calculate_size()
        actual_file.calculate_sha1()
        actual_file.has_raw = True
        actual_file.source_url = original_file.source_url
        actual_file.source_filename = original_file.source_filename
        actual_file.save()

        for sample in samples:
            original_file_sample_association = OriginalFileSampleAssociation()
            original_file_sample_association.sample = sample
            original_file_sample_association.original_file = actual_file
            original_file_sample_association.save()

        return actual_file


def download_geo(job_id: int) -> None:
    """The main function for the GEO Downloader.

    Downloads a single tar file containing the files representing
    samples relating to a single experiement stored in
    GEO.
    """
    job = utils.start_job(job_id)

    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)

    original_file = file_assocs[0].original_file
    url = original_file.source_url
    accession_code = job.accession_code

    sample_assocs = OriginalFileSampleAssociation.objects.filter(original_file=original_file)
    related_samples = Sample.objects.filter(id__in=sample_assocs.values('sample_id'))

    # First, download the sample archive URL.
    # Then, unpack all the ones downloaded.
    # Then create processor jobs!

    # The files for all of the samples are
    # contained within the same zip file. Therefore only
    # download the one.
    os.makedirs(LOCAL_ROOT_DIR + '/' + accession_code, exist_ok=True)
    dl_file_path = LOCAL_ROOT_DIR + '/' + accession_code + '/' + url.split('/')[-1]

    logger.debug("Starting to download: " + url, job_id=job_id, accession_code=accession_code)
    _download_file(url, dl_file_path, job)
    original_file.absolute_file_path = dl_file_path
    original_file.is_downloaded = True
    original_file.save()

    has_raw = True
    unpacked_sample_files = []

    # These files are tarred, and also subsequently gzipped
    if '.tar' in dl_file_path:
        try:
            extracted_files = _extract_tar(dl_file_path, accession_code)
        except Exception as e:
            job.failure_reason = e
            logger.exception(
                "Error occured while extracting tar file.", path=dl_file_path, exception=str(e))
            utils.end_downloader_job(job, success=False)
            return

        for og_file in extracted_files:

            filename = og_file['filename']
            if '_' in filename:
                sample_id = filename.split('_')[0]
            else:
                sample_id = filename.split('.')[0]

            try:
                sample = Sample.objects.get(accession_code=sample_id)
            except Exception as e:
                # We don't have this sample, but it's not a total failure. This happens.
                continue

            try:
                # Files from the GEO supplemental file are gzipped inside of the tarball. Great!
                archive_file = OriginalFile.objects.get(source_filename__contains=sample_id)
                archive_file.is_downloaded = True
                archive_file.is_archive = True
                archive_file.absolute_file_path = og_file['absolute_path']
                archive_file.calculate_size()
                archive_file.calculate_sha1()
                archive_file.save()

                if '.gz' in og_file['filename']:
                    extracted_subfile = _extract_gz(og_file['absolute_path'], accession_code)
                else:
                    extracted_subfile = [og_file]

                # Check if the OriginalFile for the file contained
                # within the archive exists already, create it if it
                # doesn't, and then check if we actually need to queue
                # it for processing or not.
                actual_file = _get_actual_file_if_queueable(extracted_subfile[0], original_file, [sample])
                if actual_file:
                    unpacked_sample_files.append(actual_file)

                archive_file.delete_local_file()
                archive_file.is_downloaded = False
                archive_file.save()
            except Exception as e:
                # TODO - is this worth failing a job for?
                logger.debug("Found a file we didn't have an OriginalFile for! Why did this happen?: "
                             + og_file['filename'],
                             exc_info=1,
                             file=og_file['filename'],
                             sample_id=sample_id,
                             accession_code=accession_code)
                # If we don't know why we have it, get rid of it.
                os.remove(og_file["absolute_path"])

    # This is a .tgz file.
    elif '.tgz' in dl_file_path:
        # If this is the MINiML file, it has been preprocessed
        if '_family.xml.tgz' in dl_file_path:
            has_raw = False

        try:
            extracted_files = _extract_tgz(dl_file_path, accession_code)
        except Exception as e:
            job.failure_reason = e
            logger.exception("Error occured while extracting tgz file.",
                             path=dl_file_path,
                             exception=str(e))
            utils.end_downloader_job(job, success=False)
            return

        for og_file in extracted_files:

            if '.txt' in og_file['filename']:
                try:
                    gsm_id = og_file['filename'].split('-')[0]
                    sample = Sample.objects.get(accession_code=gsm_id)
                except Exception as e:
                    os.remove(og_file["absolute_path"])
                    continue

                # Check if the OriginalFile for the file contained
                # within the archive exists already, create it if it
                # doesn't, and then check if we actually need to queue
                # it for processing or not.
                actual_file = _get_actual_file_if_queueable(og_file, original_file, [sample])
                if actual_file:
                    unpacked_sample_files.append(actual_file)

    # These files are only gzipped.
    # These are generally the _actually_ raw (rather than the non-raw data in a RAW file) data
    elif '.gz' in dl_file_path:
        try:
            extracted_files = _extract_gz(dl_file_path, accession_code)
        except Exception as e:
            job.failure_reason = e
            logger.exception("Error occured while extracting gz file.",
                             path=dl_file_path,
                             exception=str(e))
            utils.end_downloader_job(job, success=False)
            return

        for og_file in extracted_files:

            filename = og_file['filename']
            sample_id = filename.split('.')[0]

            try:
                # The archive we downloaded
                archive_file = OriginalFile.objects.get(source_filename__contains=filename)
                archive_file.is_downloaded = True
                archive_file.is_archive = True
                archive_file.absolute_file_path = dl_file_path
                archive_file.calculate_size()
                archive_file.calculate_sha1()
                archive_file.save()

                # Check if the OriginalFile for the file contained
                # within the archive exists already, create it if it
                # doesn't, and then check if we actually need to queue
                # it for processing or not.
                actual_file = _get_actual_file_if_queueable(og_file, original_file, related_samples)
                if actual_file:
                    unpacked_sample_files.append(actual_file)

                archive_file.delete_local_file()
                archive_file.is_downloaded = False
                archive_file.save()
            except Exception as e:
                logger.debug("Found a file we didn't have an OriginalFile for! Why did this happen?: "
                            + og_file['filename'],
                            exc_info=1,
                            file=og_file['filename'],
                            sample_id=sample_id,
                            accession_code=accession_code)
                os.remove(og_file["absolute_path"])

    # This is probably just a .txt file
    else:
        filename = dl_file_path.split('/')[-1]
        sample_id = filename.split('_')[0]

        actual_file = OriginalFile()
        actual_file.is_downloaded = True
        actual_file.is_archive = False
        actual_file.absolute_file_path = dl_file_path
        actual_file.filename = filename
        actual_file.calculate_size()
        actual_file.calculate_sha1()
        actual_file.has_raw = True
        actual_file.source_url = original_file.source_url
        actual_file.source_filename = original_file.source_filename
        actual_file.save()

        for sample in related_samples:
            new_association = OriginalFileSampleAssociation()
            new_association.original_file = actual_file
            new_association.sample = sample
            new_association.save()

        unpacked_sample_files.append(actual_file)

    if len(unpacked_sample_files) > 0:
        success = True
        logger.debug("File downloaded and extracted successfully.",
                     url=url,
                     dl_file_path=dl_file_path,
                     downloader_job=job_id)
    else:
        success = False
        logger.info("Unable to extract any files.",
                    url=url,
                    dl_file_path=dl_file_path,
                    downloader_job=job_id)
        job.failure_reason = "Failed to extract any downloaded files."

    if success:
        utils.create_processor_jobs_for_original_files(unpacked_sample_files, job)

    if original_file.is_archive:
        original_file.delete_local_file()

    utils.end_downloader_job(job, success)

    return success

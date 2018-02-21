from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import subprocess
from contextlib import closing
from data_refinery_common.models import File, DownloaderJob
from data_refinery_common.models.new_models import DownloaderJobOriginalFileAssociation
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _download_file(download_url: str, downloader_job: DownloaderJob, target_file_path: str, force_ftp: bool=False) -> bool:
    """ Download file dispatcher. Dispatches to the FTP or Aspera downloader """

    # SRA files have Apsera downloads.
    if 'ftp.sra.ebi.ac.uk' in download_url and not force_ftp:
        # From: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz
        # To: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz
        download_url = download_url.replace('ftp://', 'era-fasp@')
        download_url = download_url.replace('ftp', 'fasp')
        download_url = download_url.replace('.uk/', '.uk:/')

        return _download_file_aspera(download_url, downloader_job, target_file_path)
    else:
        return _download_file_ftp(download_url, downloader_job, target_file_path)

def _download_file_ftp(download_url: str, downloader_job: DownloaderJob, target_file_path: str) -> bool:
    """ Download a file to a location using FTP via urllib. """
    try:
        logger.debug("Downloading file from %s to %s via FTP.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)
        with closing(urllib.request.urlopen(download_url)) as request:
            with open(target_file_path, "wb") as target_file:
                shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()
    except Exception:
        logger.exception("Exception caught while downloading batch from the URL via FTP: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "batch from the URL via FTP: {}").format(download_url)
        return False

    return True

def _download_file_aspera(download_url: str, downloader_job: DownloaderJob, target_file_path: str) -> bool:
    """ Download a file to a location using Aspera by shelling out to the `ascp` client. """

    try:
        logger.debug("Downloading file from %s to %s via Aspera.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        # aspera.sra.ebi.ac.uk users port 33001 for SSH communication
        # We are also NOT using encryption (-T) to avoid slowdown,
        # and we are not using any kind of rate limiting.
        command_str = (".aspera/cli/bin/ascp -P33001 -i .aspera/cli/etc/asperaweb_id_dsa.openssh {src} {dest}")
        formatted_command = command_str.format(src=download_url,
                                               dest=target_file_path)
        completed_command = subprocess.run(formatted_command.split(),
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

        # Something went wrong! Else, just fall through to returning True.
        if completed_command.returncode != 0:
            import pdb
            pdb.set_trace()
            stderr = str(completed_command.stderr)
            logger.error("Shell call to ascp failed with error message: %s\nCommand was: %s",
                         stderr,
                         formatted_command
                        )
            return False

    except Exception:
        logger.exception("Exception caught while downloading batch from the URL via Aspera: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "batch from the URL via Aspera: {}").format(download_url)
        return False
    return True

def download_sra(job_id: int) -> None:
    """The main function for the SRA Downloader.

    Fairly straightforward, just downloads the Batch's file from SRA
    and pushes it into Temporary Storage.
    """
    job = utils.start_job(job_id)
    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)

    downloaded_files = []
    success = False
    for assoc in file_assocs:
        original_file = assoc.original_file

        if original_file.is_downloaded:
            logger.error("File already downloaded!")
            continue

        os.makedirs(LOCAL_ROOT_DIR + '/' + job.accession_code, exist_ok=True)
        dl_file_path = LOCAL_ROOT_DIR + '/' + job.accession_code + '/' + original_file.source_filename
        success = _download_file(original_file.source_url, job, dl_file_path)

        if success:
            original_file.is_downloaded = True
            original_file.absolute_file_path = dl_file_path
            original_file.file_name = original_file.source_filename
            original_file.is_archive = False
            original_file.calculate_size()
            original_file.calculate_sha1()
            original_file.save()

            downloaded_files.append(original_file)
        else:
            logger.error("A problem occured while downloading")

    utils.end_downloader_job(job, success)
    if success:
        utils.create_processor_job_for_original_files(downloaded_files)

    #         target_file_path = file.get_temp_pre_path(job_dir)
    #         success = _download_file(file, job, target_file_path)

    # batches = job.batches.all()
    # success = True
    # job_dir = utils.JOB_DIR_PREFIX + str(job_id)

    # # There should only be one batch per SRA job.
    # if batches.count() == 1:
    #     files = File.objects.filter(batch=batches[0])
    #     # All the files will be downloaded to the same directory
    #     target_directory = files[0].get_temp_dir(job_dir)
    #     os.makedirs(target_directory, exist_ok=True)
    # elif batches.count() > 1:
    #     message = "More than one batch found for SRA downloader job. There should only be one."
    #     logger.error(message, downloader_job=job_id)
    #     job.failure_reason = message
    #     success = False
    # else:
    #     message = "No batches found."
    #     logger.error(message, downloader_job=job_id)
    #     job.failure_reason = message
    #     success = False

    # if success:
    #     for file in files:
    #         target_file_path = file.get_temp_pre_path(job_dir)
    #         success = _download_file(file, job, target_file_path)

    #         # If a download fails stop the job and fail gracefully.
    #         if not success:
    #             break

    #         try:
    #             file.size_in_bytes = os.path.getsize(target_file_path)
    #             file.save()
    #             file.upload_raw_file(job_dir)
    #         except Exception:
    #             logger.exception("Exception caught while uploading file.",
    #                              downloader_job=job.id,
    #                              batch=batches[0].id,
    #                              file=file.id,
    #                              file_name=file.name)
    #             job.failure_reason = "Exception caught while uploading file."
    #             success = False
    #             break

    # if success:
    #     logger.debug("Files for batch %s downloaded and extracted successfully.",
    #                  file.download_url,
    #                  downloader_job=job_id)

    # utils.end_job(job, batches, success)

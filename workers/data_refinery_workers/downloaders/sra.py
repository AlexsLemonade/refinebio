from contextlib import closing
from typing import List
from django.utils import timezone
import os
import shutil
import subprocess
import time
import urllib.request

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
)
from data_refinery_common.utils import get_env_variable, get_fasp_sra_download
from data_refinery_workers.downloaders import utils

logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256

def _download_file(download_url: str,
                   downloader_job: DownloaderJob,
                   target_file_path: str,
                   force_ftp: bool=False) -> bool:
    """ Download file dispatcher. Dispatches to the FTP or Aspera downloader """

    # SRA files have Apsera downloads.
    if 'ftp.sra.ebi.ac.uk' in download_url and not force_ftp:
        # From: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz
        # To: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz
        download_url = download_url.replace('ftp://', 'era-fasp@')
        download_url = download_url.replace('ftp', 'fasp')
        download_url = download_url.replace('.uk/', '.uk:/')
        return _download_file_aspera(download_url, downloader_job, target_file_path, source="ENA")
    elif "ncbi.nlm.nih.gov" in download_url and not force_ftp:
        # Try to convert old-style endpoints into new-style endpoints if possible
        try:
            if 'anonftp' in download_url:
                accession = download_url.split('/')[-1].split('.sra')[0]
                new_url = get_fasp_sra_download(accession)
                if new_url:
                    download_url = new_url
        except Exception:
            pass
        return _download_file_aspera(download_url, downloader_job, target_file_path, source="NCBI")
    else:
        return _download_file_ftp(download_url, downloader_job, target_file_path)


def _download_file_ftp(download_url: str, downloader_job: DownloaderJob, target_file_path: str) -> bool:
    """ Download a file to a location using FTP via urllib. """
    try:
        logger.debug("Downloading file from %s to %s via FTP.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()

        with closing(urllib.request.urlopen(download_url)) as request:
            with open(target_file_path, "wb") as target_file:
                shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        urllib.request.urlcleanup()
    except Exception:
        logger.exception("Exception caught while downloading file from the URL via FTP: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "file from the URL via FTP: {}").format(download_url)
        return False

    return True


def _download_file_aspera(download_url: str,
                          downloader_job: DownloaderJob,
                          target_file_path: str,
                          attempt: int=0,
                          source="NCBI"
                          ) -> bool:
    """ Download a file to a location using Aspera by shelling out to the `ascp` client. """

    try:
        logger.debug("Downloading file from %s to %s via Aspera.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        if source is "ENA":
            # aspera.sra.ebi.ac.uk users port 33001 for SSH communication
            # We are also NOT using encryption (-T) to avoid slowdown,
            # and we are not using any kind of rate limiting.
            command_str = ".aspera/cli/bin/ascp -P33001 -i .aspera/cli/etc/asperaweb_id_dsa.openssh {src} {dest}"
            formatted_command = command_str.format(src=download_url,
                                                   dest=target_file_path)
            completed_command = subprocess.run(formatted_command.split(),
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        else:
            # NCBI requires encryption and recommends -k1 resume, as well as the 450m limit and -Q (play fair).
            # ex: https://github.com/AlexsLemonade/refinebio/pull/1189#issuecomment-478018580
            command_str = ".aspera/cli/bin/ascp -p -Q -T -k1 -l 450m -i .aspera/cli/etc/asperaweb_id_dsa.openssh {src} {dest}"
            formatted_command = command_str.format(src=download_url,
                                                   dest=target_file_path)
            logger.info("Starting NCBI ascp", time=str(timezone.now()))
            completed_command = subprocess.run(formatted_command.split(),
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
            logger.info("Ending NCBI ascp", time=str(timezone.now()))

        # Something went wrong! Else, just fall through to returning True.
        if completed_command.returncode != 0:

            stdout = completed_command.stdout.decode().strip()
            stderr = completed_command.stderr.decode().strip()
            logger.debug("Shell call of `%s` to ascp failed with error message: %s",
                         formatted_command,
                         stderr,
                         downloader_job=downloader_job.id)

            # Sometimes, Aspera fails mysteriously.
            # Wait a few minutes and try again.
            if attempt > 5:
                logger.info("Final shell call of `%s` to ascp failed with error message: %s",
                         formatted_command,
                         stderr + "\nSTDOUT: " + stdout,
                         downloader_job=downloader_job.id)
                downloader_job.failure_reason = "stderr:\n " + stderr + "\nstdout:\n " + stdout
                return False
            else:
                time.sleep(5)
                return _download_file_aspera(download_url,
                                             downloader_job,
                                             target_file_path,
                                             attempt + 1,
                                             source
                                             )
    except Exception:
        logger.exception("Exception caught while downloading file from the URL via Aspera: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "file from the URL via Aspera: {}").format(download_url)
        return False

    # If Aspera has given a zero-byte file for some reason, let's back off and retry.
    if (not os.path.exists(target_file_path)) or (os.path.getsize(target_file_path) < 1):
        if os.path.exists(target_file_path):
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
                                     attempt + 1,
                                     source
                                     )
    return True


def download_sra(job_id: int) -> None:
    """The main function for the SRA Downloader.

    Fairly straightforward, just downloads the file from SRA.
    """
    job = utils.start_job(job_id)
    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)

    downloaded_files = []
    success = None
    for assoc in file_assocs:
        original_file = assoc.original_file

        if original_file.is_downloaded:
            logger.info("File already downloaded!",
                         original_file_id=original_file.id,
                         downloader_job=job_id)
            success = True
            continue

        sample_accession_code = original_file.samples.first().accession_code
        exp_path = LOCAL_ROOT_DIR + '/' + job.accession_code
        samp_path = exp_path + '/' + sample_accession_code
        os.makedirs(exp_path, exist_ok=True)
        os.makedirs(samp_path, exist_ok=True)
        dl_file_path = samp_path + '/' + original_file.source_filename
        success = _download_file(original_file.source_url, job, dl_file_path)

        if success:
            original_file.is_downloaded = True
            original_file.absolute_file_path = dl_file_path
            original_file.filename = original_file.source_filename
            original_file.is_archive = False
            original_file.calculate_size()
            original_file.calculate_sha1()
            original_file.save()

            downloaded_files.append(original_file)
        else:
            break

    if success:
        utils.create_processor_job_for_original_files(downloaded_files, job)

    utils.end_downloader_job(job, success)

    return success, downloaded_files

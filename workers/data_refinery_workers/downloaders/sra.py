from contextlib import closing
from typing import List
from django.utils import timezone
from ftplib import FTP
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
    OriginalFileSampleAssociation,
    Sample,
)
from data_refinery_common.rna_seq import _build_ena_file_url
from data_refinery_common.utils import get_env_variable, get_https_sra_download
from data_refinery_workers.downloaders import utils
from data_refinery_common.job_management import create_processor_job_for_original_files

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
            if 'anonftp' in download_url or 'dbtest' in download_url:
                accession = download_url.split('/')[-1].split('.sra')[0]
                new_url = get_https_sra_download(accession)
                if new_url:
                    download_url = new_url
        except Exception:
            pass
        return _download_file_http(download_url, downloader_job, target_file_path)
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


def _download_file_http(download_url: str,
                        downloader_job: DownloaderJob,
                        target_file_path: str
                          ) -> bool:
    try:
        target_file = open(target_file_path, "wb")
        logger.debug("Downloading file from %s to %s using HTTP.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        with closing(urllib.request.urlopen(download_url, timeout=60)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logger.exception("Exception caught while downloading file.",
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = "Exception caught while downloading file" \
                                        + str(e).replace('\n', '\\n')
        return False
    finally:
        target_file.close()

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


def _has_unmated_reads(accession_code: str) -> bool:
    """Checks if the SRA accession has unmated reads.

    Returns True if it does and False if it doesn't."""
    full_ftp_link = _build_ena_file_url(accession_code)

    # Strip off the protocol code because we know it's FTP and the FTP
    # library doesn't want the protocol.
    no_protocol_link = full_ftp_link.split("://")[1]

    # We need to extract the server so we can login to it.
    split_link = no_protocol_link.split("/")
    ftp_server = split_link[0]

    # We need to get the FTP directory the file is in so we can check
    # how many other files are in it. Therefore we're looking to get
    # the path between the server and the filename itself.
    sample_directory = "/".join(split_link[1:-1])

    ftp = None
    try:
        ftp = FTP(ftp_server)
        ftp.login()
        ftp.cwd(sample_directory)

        # If there's three files then there's unmated reads, because
        # there's the one read file, the other read file, and the
        # unmated reads.
        if len(ftp.nlst()) == 3:
            return True
        else:
            return False
    except:
        # If we can't find the sample on ENA's FTP server, then we
        # shouldn't try to download it from there.
        return False
    finally:
        if ftp:
            ftp.close()

    # Shouldn't reach here, but just in case default to NCBI.
    return False


def _replace_dotsra_with_fastq_files(sample: Sample,
                                     downloader_job: DownloaderJob,
                                     original_file: OriginalFile) -> List[OriginalFile]:
    """Replaces a .SRA file with two .fastq files.

    This function should only be called on a sample which has unmated
    reads, so it makes the assumption that the sample passed into it
    has at least two read files in ENA.
    """
    read_one_url = _build_ena_file_url(sample.accession_code, "_1")
    read_two_url = _build_ena_file_url(sample.accession_code, "_2")

    # Technically this is a different file, but deleting this one and
    # its associations just to recreate another with the same
    # associations seems rather pointless.
    original_file.source_url = read_one_url
    original_file.source_filename = read_one_url.split('/')[-1]
    original_file.save()

    read_two_original_file = OriginalFile.objects.get_or_create(
        source_url = read_two_url,
        source_filename = read_two_url.split('/')[-1],
        has_raw = True
    )[0]
    OriginalFileSampleAssociation.objects.get_or_create(
        original_file = read_two_original_file,
        sample = sample
    )
    DownloaderJobOriginalFileAssociation.objects.get_or_create(
        original_file = read_two_original_file,
        downloader_job = downloader_job
    )
    return [original_file, read_two_original_file]


def download_sra(job_id: int) -> None:
    """The main function for the SRA Downloader.

    Fairly straightforward, just downloads the file from SRA.
    """
    job = utils.start_job(job_id)
    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
    original_files = job.original_files.all()

    original_file = original_files[0]
    sample = original_file.samples.first()
    if _has_unmated_reads(sample.accession_code):
        original_files = _replace_dotsra_with_fastq_files(sample, job, original_file)
    else:
        # _replace_dotsra_with_fastq_files returns a list of
        # OriginalFiles so turn the queryset of
        # DownloaderJobOriginalFileAssociations into a list of
        # OriginalFiles to be consistent.
        original_files = [assoc.original_file for assoc in file_assocs]

    downloaded_files = []
    success = None
    for original_file in original_files:
        if original_file.is_downloaded:
            logger.info("File already downloaded!",
                         original_file_id=original_file.id,
                         downloader_job=job_id)
            success = True
            continue

        exp_path = LOCAL_ROOT_DIR + '/' + job.accession_code
        samp_path = exp_path + '/' + sample.accession_code
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
        create_processor_job_for_original_files(downloaded_files, job)

    utils.end_downloader_job(job, success)

    return success, downloaded_files

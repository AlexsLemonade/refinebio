from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
from contextlib import closing
from data_refinery_common.models import File, DownloaderJob
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _download_file(file: File, downloader_job: DownloaderJob, target_file_path: str) -> bool:
    try:
        logger.debug("Downloading file from %s to %s.",
                     file.download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)
        with closing(urllib.request.urlopen(file.download_url)) as request:
            with open(target_file_path, "wb") as target_file:
                shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()
    except Exception:
        logger.exception("Exception caught while downloading batch from the URL: %s",
                         file.download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "batch from the URL: {}").format(file.download_url)
        return False

    return True


def download_sra(job_id: int) -> None:
    job = utils.start_job(job_id)
    batches = job.batches.all()
    success = True
    job_dir = utils.JOB_DIR_PREFIX + str(job_id)

    # There should only be one batch per SRA job.
    if batches.count() == 1:
        files = File.objects.filter(batch=batches[0])
        # All the files will be downloaded to the same directory
        target_directory = files[0].get_temp_dir(job_dir)
        os.makedirs(target_directory, exist_ok=True)
    elif batches.count() > 1:
        message = "More than one batch found for SRA downloader job. There should only be one."
        logger.error(message, downloader_job=job_id)
        job.failure_reason = message
        success = False
    else:
        message = "No batches found."
        logger.error(message, downloader_job=job_id)
        job.failure_reason = message
        success = False

    if success:
        for file in files:
            target_file_path = file.get_temp_pre_path(job_dir)
            success = _download_file(file, job, target_file_path)

            # If a download fails stop the job and fail gracefully.
            if not success:
                break

            try:
                file.size_in_bytes = os.path.getsize(target_file_path)
                file.save()
                file.upload_raw_file(job_dir)
            except Exception:
                logger.exception("Exception caught while uploading file.",
                                 downloader_job=job.id,
                                 batch=batches[0].id,
                                 file=file.id,
                                 file_name=file.name)
                job.failure_reason = "Exception caught while uploading file."
                success = False
                break

    if success:
        logger.debug("Files for batch %s downloaded and extracted successfully.",
                     file.download_url,
                     downloader_job=job_id)

    utils.end_job(job, batches, success)

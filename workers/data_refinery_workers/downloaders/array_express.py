from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import zipfile
from typing import List
from contextlib import closing
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_common.models import Batch, DownloaderJob
from data_refinery_common import file_management
from data_refinery_workers.downloaders import utils
import logging


logger = get_task_logger(__name__)


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256
JOB_DIR_PREFIX = "downloader_job_"


def _verify_batch_grouping(batches: List[Batch], job: DownloaderJob) -> None:
    """All batches in the same job should have the same downloader url"""
    for batch in batches:
        if batch.download_url != batches[0].download_url:
            failure_message = "A Batch doesn't have the same download URL as the other batches"
            logger.error(failure_message + " in Downloader Job #%d.",
                         job.id)
            job.failure_reason = failure_message
            raise ValueError(failure_message)


def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    try:
        logger.debug("Downloading file from %s to %s. (Job #%d)",
                     download_url,
                     file_path,
                     job.id)
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logging.exception("Exception caught while downloading batch in Downloader Job #%d.",
                          job.id)
        job.failure_reason = "Exception caught while downloading batch"
        raise
    finally:
        target_file.close()


def _extract_file(batches: List[Batch], job: DownloaderJob) -> None:
    """Extract zip from temp directory and move to raw directory.

    Additionally this function sets the size_in_bytes field of each
    Batch in batches. To save database calls it does not save the
    batch itself since it will be saved soon when its status
    changes in utils.end_job.
    """
    # zip_path and local_dir should be common to all batches in the group
    job_dir = JOB_DIR_PREFIX + str(job.id)
    zip_path = file_management.get_temp_download_path(batches[0], job_dir)
    local_dir = file_management.get_temp_dir(batches[0], job_dir)
    dirs_to_clean = set()

    logger.debug("Extracting %s for Downloader Job %d.", zip_path, job.id)

    try:
        zip_ref = zipfile.ZipFile(zip_path, "r")
        zip_ref.extractall(local_dir)

        for batch in batches:
            batch_directory = file_management.get_temp_dir(batch, job_dir)
            raw_file_location = file_management.get_temp_pre_path(batch, job_dir)

            # The platform is part of the batch's location so if the
            # batches in this job have different platforms then some
            # of them need to be moved to the directory corresponding
            # to thier platform.
            if local_dir != batch_directory:
                os.makedirs(batch_directory, exist_ok=True)
                dirs_to_clean.add(batch_directory)
                incorrect_location = os.path.join(local_dir, batch.name)
                os.rename(incorrect_location, raw_file_location)

            batch.size_in_bytes = os.path.getsize(raw_file_location)
            file_management.upload_raw_file(batch, job_dir)
    except Exception:
        logging.exception("Exception caught while extracting %s during Downloader Job #%d.",
                          zip_path,
                          job.id)
        job.failure_reason = "Exception caught while extracting " + zip_path
        raise
    finally:
        zip_ref.close()
        file_management.remove_temp_directory(batches[0], job_dir)
        for directory in dirs_to_clean:
            shutil.rmtree(directory)


@shared_task
def download_array_express(job_id: int) -> None:
    job = utils.start_job(job_id)
    batches = job.batches.all()
    success = True
    job_dir = JOB_DIR_PREFIX + str(job_id)

    if batches.count() > 0:
        target_directory = file_management.get_temp_dir(batches[0], job_dir)
        os.makedirs(target_directory, exist_ok=True)
        target_file_path = file_management.get_temp_download_path(batches[0], job_dir)
        download_url = batches[0].download_url
    else:
        logger.error("No batches found for Downloader Job #%d.",
                     job_id)
        success = False

    if success:
        try:
            _verify_batch_grouping(batches, job)

            # The files for all of the batches in the grouping are
            # contained within the same zip file. Therefore only
            # download the one.
            _download_file(download_url, target_file_path, job)
            _extract_file(batches, job)
        except Exception:
            # Exceptions are already logged and handled.
            # Just need to mark the job as failed.
            success = False

    if success:
        logger.debug("File %s downloaded and extracted successfully in Downloader Job #%d.",
                     download_url,
                     job_id)

    utils.end_job(job, batches, success)

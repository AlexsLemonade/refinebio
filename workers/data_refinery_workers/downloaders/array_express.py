from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import zipfile
from typing import List
from contextlib import closing
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_models.models import Batch
from data_refinery_common import file_management
from data_refinery_workers.downloaders import utils
import logging


logger = get_task_logger(__name__)


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256
JOB_DIR_PREFIX = "downloader_job_"


def _verify_batch_grouping(batches: List[Batch], job_id: int) -> None:
    """All batches in the same job should have the same downloader url"""
    for batch in batches:
        if batch.download_url != batches[0].download_url:
            logger.error(("A Batch doesn't have the same download URL as the other batches"
                          " in downloader job #%d."),
                         job_id)
            raise ValueError("A batch doesn't have the same download url as other batches.")


def _download_file(download_url: str, file_path: str, job_id: int) -> None:
    try:
        logger.debug("Downloading file from %s to %s. (Job #%d)",
                     download_url,
                     file_path,
                     job_id)
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logging.exception("Exception caught while running Job #%d.",
                          job_id)
        raise
    finally:
        target_file.close()


def _extract_file(batches: List[Batch], job_id: int) -> None:
    """Extract zip from temp directory and move to raw directory."""
    # zip_path and local_dir should be common to all batches in the group
    job_dir = JOB_DIR_PREFIX + str(job_id)
    zip_path = file_management.get_temp_download_path(batches[0], job_dir)
    local_dir = file_management.get_temp_dir(batches[0], job_dir)

    logger.debug("Extracting %s for job %d.", zip_path, job_id)

    try:
        zip_ref = zipfile.ZipFile(zip_path, "r")
        zip_ref.extractall(local_dir)

        for batch in batches:
            file_management.upload_raw_file(batch, job_dir)
    except Exception:
        logging.exception("Exception caught while extracting %s during Job #%d.",
                          zip_path,
                          job_id)
        raise
    finally:
        zip_ref.close()
        file_management.remove_temp_directory(batches[0], job_dir)


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
        logger.error("No batches found for job #%d.",
                     job_id)
        success = False

    if success:
        try:
            _verify_batch_grouping(batches, job_id)

            # The files for all of the batches in the grouping are
            # contained within the same zip file. Therefore only
            # download the one.
            _download_file(download_url, target_file_path, job_id)
            _extract_file(batches, job_id)
        except Exception:
            # Exceptions are already logged and handled.
            # Just need to mark the job as failed.
            success = False

    if success:
        logger.debug("File %s downloaded and extracted successfully in Job #%d.",
                     download_url,
                     job_id)

    utils.end_job(job, batches, success)

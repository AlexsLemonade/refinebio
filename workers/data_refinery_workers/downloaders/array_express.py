from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import zipfile
from typing import List
from contextlib import closing
from celery import shared_task
from celery.utils.log import get_task_logger
from django.core.exceptions import ObjectDoesNotExist
from data_refinery_models.models import (
    Batch,
    DownloaderJob,
    DownloaderJobsToBatches
)
from data_refinery_workers.downloaders import utils
import logging

logger = get_task_logger(__name__)

# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _verify_batch_grouping(batches: List[Batch], job_id):
    """All batches in the same job should have the same downloader url"""
    for batch in batches:
        if batch.download_url != batches[0].download_url:
            logger.error(("A Batch doesn't have the same download URL as the other batches"
                          " in downloader job #%d."),
                         job_id)
            raise ValueError("A batch doesn't have the same download url as other batches.")


def _download_file(download_url, file_path, job_id):
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


def _extract_file(file_path, job_id):
    try:
        zip_ref = zipfile.ZipFile(file_path, 'r')
        zip_ref.extractall(os.path.dirname(file_path))
    except Exception:
        logging.exception("Exception caught while extracting %s during Job #%d.",
                          file_path,
                          job_id)
        raise
    finally:
        zip_ref.close()
        os.remove(file_path)


@shared_task
def download_array_express(job_id):
    logger.debug("Starting job with id: %s.", job_id)
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except ObjectDoesNotExist:
        logger.error("Cannot find downloader job record with ID %d.", job_id)
        return

    success = True
    utils.start_job(job)

    batch_relations = DownloaderJobsToBatches.objects.filter(downloader_job_id=job_id)
    batches = [br.batch for br in batch_relations]

    if len(batches) > 0:
        target_file_path = utils.prepare_destination(batches[0])
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
            _extract_file(target_file_path, job_id)
        except Exception:
            # Exceptions are already logged and handled.
            # Just need to mark the job as failed.
            success = False

    if success:
        logger.debug("File %s downloaded and extracted successfully in Job #%d.",
                     download_url,
                     job_id)

    utils.end_job(job, batches, success)

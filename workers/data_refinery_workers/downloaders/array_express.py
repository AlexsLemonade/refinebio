from __future__ import absolute_import, unicode_literals
import requests
from celery import shared_task
from celery.utils.log import get_task_logger
from django.core.exceptions import ObjectDoesNotExist
from data_refinery_models.models import Batch, DownloaderJob
from data_refinery_workers.downloaders import utils

logger = get_task_logger(__name__)

# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


@shared_task
def download_array_express(job_id):
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except ObjectDoesNotExist:
        logger.error("Cannot find downloader job record with ID %d.", job_id)
        return

    success = True
    utils.start_job(job)

    try:
        batch = Batch.objects.get(id=job.batch_id)
    except ObjectDoesNotExist:
        logger.error("Cannot find batch record with ID %d.", job.batch_id)
        utils.end_job(job, None, False)
        return

    target_file_path = utils.prepare_destination(batch)

    logger.info("Downloading file from %s to %s. (Batch #%d, Job #%d)",
                batch.download_url,
                target_file_path,
                batch.id,
                job_id)

    try:
        target_file = open(target_file_path, "wb")
        request = requests.get(batch.download_url, stream=True)

        for chunk in request.iter_content(CHUNK_SIZE):
            if chunk:
                target_file.write(chunk)
                target_file.flush()
    except Exception as e:
        success = False
        logger.error("Exception caught while running Job #%d for Batch #%d "
                     + "with message: %s",
                     job_id,
                     batch.id,
                     e)
    finally:
        target_file.close()
        request.close()

    if success:
        logger.info("File %s (Batch #%d) downloaded successfully in Job #%d.",
                    batch.download_url,
                    batch.id,
                    job_id)

    utils.end_job(job, batch, success)

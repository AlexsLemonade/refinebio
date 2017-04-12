from __future__ import absolute_import, unicode_literals
import os
import requests
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_models.models import Batch, DownloaderJob
from data_refinery_workers.downloaders import utils

logger = get_task_logger(__name__)

# chunk_size is in bytes
CHUNK_SIZE = 4096
# This path is within the Docker container.
ROOT_URI = "/home/user/data_store/"


@shared_task
def download_array_express(job_id):
    job = (DownloaderJob
           .objects
           .filter(id=job_id)
           [:1]
           .get())

    utils.start_job(job)

    batch = (Batch
             .objects
             .filter(id=job.batch_id)
             [:1]
             .get())

    success = True
    target_directory = ROOT_URI + batch.internal_location
    os.makedirs(target_directory, exist_ok=True)

    filename = batch.download_url.split('/')[-1]
    target_file_name = target_directory + filename

    logger.info("Downloading file from %s to %s. (Batch #%d, Job #%d)",
                batch.download_url,
                target_file_name,
                batch.id,
                job_id)

    try:
        target_file = open(target_file_name, "wb")
        request = requests.get(batch.download_url, stream=True)

        raise Exception("test")
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

    if success:
        logger.info("File %s (Batch #%d) downloaded successfully in Job #%d.",
                    batch.download_url,
                    batch.id,
                    job_id)

    utils.end_job(job, batch, success)

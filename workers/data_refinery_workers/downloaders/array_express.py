from __future__ import absolute_import, unicode_literals
import os
import requests
from celery import shared_task
from celery.utils.log import get_task_logger
from django.utils import timezone
from data_refinery_models.models import Batch, DownloaderJob

logger = get_task_logger(__name__)

# chunk_size is in bytes
CHUNK_SIZE = 4096
# This path is within the Docker container
ROOT_URI = "/home/user/data_store/"


def start_job(job: DownloaderJob):
    """Record in the database that this job is being started.
    This should be moved to a more reusable location."""
    job.worker_id = "For now there's only one. For now..."
    job.start_time = timezone.now()
    job.save()


def end_job(job: DownloaderJob, batch: Batch):
    """Record in the database that this job has completed successfully.
    This should be moved to a more reusable location.
    This should also queue a processor job for the batch.
    Should it? That may not be entirely doable from batch info alone..."""
    job.success = True
    job.end_time = timezone.now()
    job.save()


@shared_task
def download_array_express(job_id):
    job = (DownloaderJob
           .objects
           .filter(id=job_id)
           [:1]
           .get())

    start_job(job)

    batch = (Batch
             .objects
             .filter(id=job.batch_id)
             [:1]
             .get())

    target_directory = ROOT_URI + batch.internal_location
    os.makedirs(target_directory, exist_ok=True)

    filename = batch.download_url.split('/')[-1]
    target_file_name = target_directory + filename

    logger.info("Downloading file from %s to %s. (Batch #%d)",
                batch.download_url,
                target_file_name,
                batch.id)

    target_file = open(target_file_name, "wb")
    request = requests.get(batch.download_url, stream=True)

    for chunk in request.iter_content(CHUNK_SIZE):
        if chunk:
            target_file.write(chunk)
            target_file.flush()

    target_file.close()
    end_job(job, batch)

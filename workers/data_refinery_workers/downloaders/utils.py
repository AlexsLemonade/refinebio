import os
import urllib
from retrying import retry
from django.utils import timezone
from data_refinery_models.models import (
    Batch,
    BatchStatuses,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_workers.processors.processor_registry \
    import processor_pipeline_registry

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# This path is within the Docker container.
ROOT_URI = "/home/user/data_store/raw"

# This path is within the Docker container.
ROOT_URI = "/home/user/data_store/raw"


def start_job(job: DownloaderJob):
    """Record in the database that this job is being started. """
    job.worker_id = "For now there's only one. For now..."
    job.start_time = timezone.now()
    job.save()


def end_job(job: DownloaderJob, batch: Batch, success):
    """Record in the database that this job has completed,
    create a processor job, and queue a processor task."""
    job.success = success
    job.end_time = timezone.now()
    job.save()

    @retry(stop_max_attempt_number=3)
    def save_batch_create_job():
        batch.status = BatchStatuses.DOWNLOADED.value
        batch.save()

        logger.info("Creating processor job for batch #%d.", batch.id)
        processor_job = ProcessorJob(batch=batch)
        processor_job.save()
        return processor_job

    @retry(stop_max_attempt_number=3)
    def queue_task(processor_job):
        processor_task = processor_pipeline_registry[batch.pipeline_required]
        processor_task.delay(processor_job.id)

    if batch is not None:
        processor_job = save_batch_create_job()
        queue_task(processor_job)


def prepare_destination(batch: Batch):
    """Prepare the destination directory and return the full
    path the Batch's file should be downloaded to."""
    target_directory = os.path.join(ROOT_URI, batch.internal_location)
    os.makedirs(target_directory, exist_ok=True)

    path = urllib.parse.urlparse(batch.download_url).path
    filename = os.path.basename(path)
    return os.path.join(target_directory, filename)

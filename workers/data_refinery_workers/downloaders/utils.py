import os
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
ROOT_URI = "/home/user/data_store/raw/"


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

    if batch is not None:
        batch.status = BatchStatuses.DOWNLOADED.value
        batch.save()

    logger.info("Creating processor job for batch #%d.", batch.id)
    processor_job = ProcessorJob(batch=batch)
    processor_job.save()
    processor_task = processor_pipeline_registry[batch.pipeline_required]
    processor_task.delay(processor_job.id)


def prepare_destination(batch: Batch):
    """Prepare the destination directory and return the full
    path the Batch's file should be downloaded to."""
    target_directory = ROOT_URI + batch.internal_location
    os.makedirs(target_directory, exist_ok=True)

    filename = batch.download_url.split('/')[-1]
    return target_directory + filename

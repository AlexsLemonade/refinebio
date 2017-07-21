from retrying import retry
from django.utils import timezone
from django.db import transaction
from data_refinery_common.utils import get_worker_id
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


def start_job(job_id: int) -> DownloaderJob:
    """Record in the database that this job is being started.

    Retrieves the job from the database and returns it after marking
    it as started.
    """
    logger.info("Starting job with id: %s.", job_id)
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except DownloaderJob.DoesNotExist:
        logger.error("Cannot find downloader job record with ID %d.", job_id)
        raise

    job.worker_id = get_worker_id()
    job.start_time = timezone.now()
    job.save()

    return job


def end_job(job: DownloaderJob, batches: Batch, success):
    """Record in the database that this job has completed.

    Create a processor job and queue a processor task for each batch
    if the job was successful.
    """
    @retry(stop_max_attempt_number=3)
    def save_batch_create_job(batch):
        batch.status = BatchStatuses.DOWNLOADED.value
        batch.save()

        logger.debug("Creating processor job for batch #%d.", batch.id)
        processor_job = ProcessorJob.create_job_and_relationships(
            batches=[batch], pipeline_applied=batch.pipeline_required)
        return processor_job

    @retry(stop_max_attempt_number=3)
    def queue_task(processor_job):
        processor_task = processor_pipeline_registry[batch.pipeline_required]
        processor_task.delay(processor_job.id)

    if success:
        for batch in batches:
            with transaction.atomic():
                processor_job = save_batch_create_job(batch)
                queue_task(processor_job)

    job.success = success
    job.end_time = timezone.now()
    job.save()

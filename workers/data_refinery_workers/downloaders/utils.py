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
from data_refinery_workers.task_runner import app
from data_refinery_workers._version import __version__
from data_refinery_common.job_lookup import ProcessorPipeline, PROCESSOR_PIPELINE_LOOKUP
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


def start_job(job_id: int) -> DownloaderJob:
    """Record in the database that this job is being started.

    Retrieves the job from the database and returns it after marking
    it as started.
    """
    logger.info("Starting Downloader Job.", downloader_job=job_id)
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except DownloaderJob.DoesNotExist:
        logger.error("Cannot find downloader job record.", downloader_job=job_id)
        raise

    job.worker_id = get_worker_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    return job


def end_job(job: DownloaderJob, batches: Batch, success: bool):
    """Record in the database that this job has completed.

    Create a processor job and queue a processor task for each batch
    if the job was successful.
    """
    @retry(stop_max_attempt_number=3)
    def save_batch_create_job(batch):
        batch.status = BatchStatuses.DOWNLOADED.value
        batch.save()

        # TEMPORARY for Jackie's grant:
        if batch.pipeline_required != ProcessorPipeline.NONE.value:
            logger.debug("Creating processor job for Batch.",
                         downloader_job=job.id,
                         batch=batch.id)
            with transaction.atomic():
                processor_job = ProcessorJob.create_job_and_relationships(
                    batches=[batch], pipeline_applied=batch.pipeline_required)
            return processor_job
        else:
            logger.debug("Not queuing a processor job for batch.",
                         downloader_job=job.id,
                         batch=batch.id)
            return None

    @retry(stop_max_attempt_number=3)
    def queue_task(processor_job: ProcessorJob, batch: Batch):
        if batch.pipeline_required in PROCESSOR_PIPELINE_LOOKUP:
            processor_task = PROCESSOR_PIPELINE_LOOKUP[batch.pipeline_required]
            app.send_task(processor_task, args=[processor_job.id])
            return True
        else:
            failure_template = "Could not find Processor Pipeline {} in the lookup."
            failure_message = failure_template.format(batch.pipeline_required)
            logger.error(failure_message, downloader_job=job.id, batch=batch.id)
            processor_job.failure_reason = failure_message
            processor_job.success = False
            processor_job.retried = True
            processor_job.save()
            return False

    if success:
        for batch in batches:
            processor_job = save_batch_create_job(batch)
            if batch.pipeline_required != ProcessorPipeline.NONE.value:
                try:
                    success = queue_task(processor_job, batch)
                except:
                    # If the task doesn't get sent we don't want the
                    # processor_job to be left floating
                    processor_job.delete()

                    success = False
                    job.failure_message = "Could not queue processor job task."
                    logger.error(job.failure_message)

                if success:
                    logger.info("Downloader job completed successfully.", downloader_job=job.id)

    job.success = success
    job.end_time = timezone.now()
    job.save()

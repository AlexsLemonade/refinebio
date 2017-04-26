from django.utils import timezone
from data_refinery_models.models import Batch, BatchStatuses, DownloaderJob


def start_job(job: DownloaderJob):
    """Record in the database that this job is being started. """
    job.worker_id = "For now there's only one. For now..."
    job.start_time = timezone.now()
    job.save()


def end_job(job: DownloaderJob, batch: Batch, success):
    """Record in the database that this job has completed.
    This should also queue a processor job at some point."""
    job.success = success
    job.end_time = timezone.now()
    job.save()

    if batch is not None:
        batch.status = BatchStatuses.DOWNLOADED.value
        batch.save()

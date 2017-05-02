import os
import urllib
from django.utils import timezone
from data_refinery_models.models import Batch, BatchStatuses, DownloaderJob

# This path is within the Docker container.
ROOT_URI = "/home/user/data_store/raw"


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


def prepare_destination(batch: Batch):
    """Prepare the destination directory and return the full
    path the Batch's file should be downloaded to."""
    target_directory = os.path.join(ROOT_URI, batch.internal_location)
    os.makedirs(target_directory, exist_ok=True)

    path = urllib.parse.urlparse(batch.download_url).path
    filename = os.path.basename(path)
    return os.path.join(target_directory, filename)

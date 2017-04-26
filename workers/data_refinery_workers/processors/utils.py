import os
import shutil
from django.utils import timezone
from django.core.exceptions import ObjectDoesNotExist
from typing import List, Dict, Callable
from data_refinery_models.models import Batch, BatchStatuses, ProcessorJob

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# This path is within the Docker container.
ROOT_URI = "/home/user/data_store/"


def start_job(kwargs: Dict):
    """Record in the database that this job is being started and
    retrieve the job and batch from the database.
    Retrieves the job and the job's batch from the database and
    adds them to the dictionary passed in with the keys 'job'
    and 'batch' respectively."""
    job = kwargs["job"]
    job.worker_id = "For now there's only one. For now..."
    job.start_time = timezone.now()
    job.save()

    try:
        batch = Batch.objects.get(id=job.batch_id)
    except ObjectDoesNotExist:
        logger.error("Cannot find batch record with ID %d.", job.batch_id)
        return {"success": False}

    kwargs["batch"] = batch
    return kwargs


def end_job(kwargs: Dict):
    """Record in the database that this job has completed and that
    the batch has been processed if successful."""
    job = kwargs["job"]

    if "success" in kwargs:
        success = kwargs["success"]
    else:
        success = True

    job.success = success
    job.end_time = timezone.now()
    job.save()

    if job.success:
        batch = kwargs["batch"]
        batch.status = BatchStatuses.PROCESSED.value
        batch.save()

    # Every processor returns a dict, however end_job is always called
    # last so it doesn't need to contain anything.
    return {}


def cleanup_temp_data(kwargs: Dict):
    """Removes data from raw/ and temp/ directories related to the batch."""
    batch = kwargs["batch"]

    raw_file_name = batch.download_url.split('/')[-1]
    raw_file_location = (ROOT_URI + "raw/" + batch.internal_location
                         + raw_file_name)
    temp_directory = ROOT_URI + "temp/" + batch.internal_location
    os.remove(raw_file_location)
    shutil.rmtree(temp_directory)

    return kwargs


def run_pipeline(start_value: Dict, pipeline: List[Callable]):
    """Runs a pipeline of processor functions.

    start_value must contain a key 'job_id' which is a valid id
    for a ProcessorJob record.

    Each processor fuction must accept a dictionary and return a
    dictionary.

    Any processor function which returns a dictionary
    containing a key of 'success' with a value of False will cause
    the pipeline to terminate with a call to utils.end_job.

    The key 'job' is reserved for the ProcessorJob currently being run.
    The key 'batch' is reserved for the Batch that is currently being
    processed.
    It is required that the dictionary returned by each processor
    function preserve the mappings for 'job' and 'batch' that were
    passed into it.
    """

    job_id = start_value["job_id"]
    try:
        job = ProcessorJob.objects.get(id=job_id)
    except ObjectDoesNotExist:
        logger.error("Cannot find processor job record with ID %d.", job_id)
        return

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        last_result = processor(last_result)
        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor %s failed. Terminating pipeline.",
                         processor.__name__)
            end_job(last_result)
            break

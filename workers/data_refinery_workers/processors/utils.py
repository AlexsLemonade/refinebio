from typing import List, Dict, Callable
from django.utils import timezone
from data_refinery_models.models import BatchStatuses, ProcessorJob
from data_refinery_common import file_management
from data_refinery_common.utils import get_worker_id
from data_refinery_workers._version import __version__

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def start_job(kwargs: Dict):
    """A processor function to start jobs.

    Record in the database that this job is being started and
    retrieves the job's batches from the database and adds them to the
    dictionary passed in with the key 'batches'.
    """
    job = kwargs["job"]
    job.worker_id = get_worker_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    batches = list(job.batches.all())

    if len(batches) == 0:
        logger.error("No batches found for job #%d.", job.id)
        return {"success": False}

    kwargs["batches"] = batches
    return kwargs


def end_job(kwargs: Dict):
    """A processor function to end jobs.

    Record in the database that this job has completed and that
    the batch has been processed if successful.
    """
    job = kwargs["job"]

    if "success" in kwargs:
        success = kwargs["success"]
    else:
        success = True

    job.success = success
    job.end_time = timezone.now()
    job.save()

    if job.success:
        batches = kwargs["batches"]
        for batch in batches:
            batch.status = BatchStatuses.PROCESSED.value
            batch.save()

    # Every processor returns a dict, however end_job is always called
    # last so it doesn't need to contain anything.
    return {}


def upload_processed_files(kwargs: Dict) -> Dict:
    """Uploads all the processed files for the job."""
    for batch in kwargs["batches"]:
        try:
            file_management.upload_processed_file(batch)
        except Exception:
            logging.exception(("Exception caught while uploading processed file %s for batch %d"
                               " during Job #%d."),
                              file_management.get_temp_post_path(batch),
                              batch.id,
                              kwargs["job_id"])
            processed_name = file_management.get_processed_name(batch)
            failure_template = "Exception caught while uploading processed file {}"
            kwargs["job"].failure_reason = failure_template.format(processed_name)
            kwargs["success"] = False
            return kwargs
        finally:
            file_management.remove_temp_directory(batch)

    return kwargs


def cleanup_raw_files(kwargs: Dict) -> Dict:
    """Tries to clean up raw files for the job.

    If we fail to remove the raw files, the job is still done enough
    to call a success, therefore we don't mark it as a failure.
    However logging will be important so the problem can be
    identified and the raw files cleaned up.
    """
    for batch in kwargs["batches"]:
        try:
            file_management.remove_raw_files(batch)
        except:
            # If we fail to remove the raw files, the job is still done
            # enough to call a success. However logging will be important
            # so the problem can be identified and the raw files cleaned up.
            logging.exception(("Exception caught while removing raw files %s for batch %d"
                               " during Job #%d."),
                              file_management.get_temp_pre_path(batch),
                              batch.id,
                              kwargs["job_id"])

    return kwargs


def run_pipeline(start_value: Dict, pipeline: List[Callable]):
    """Runs a pipeline of processor functions.

    start_value must contain a key 'job_id' which is a valid id for a
    ProcessorJob record.

    Each processor fuction must accept a dictionary and return a
    dictionary.

    Any processor function which returns a dictionary containing a key
    of 'success' with a value of False will cause the pipeline to
    terminate with a call to utils.end_job.

    The key 'job' is reserved for the ProcessorJob currently being
    run.  The key 'batches' is reserved for the Batches that are
    currently being processed.  It is required that the dictionary
    returned by each processor function preserve the mappings for
    'job' and 'batches' that were passed into it.
    """

    job_id = start_value["job_id"]
    try:
        job = ProcessorJob.objects.get(id=job_id)
    except ProcessorJob.DoesNotExist:
        logger.error("Cannot find processor job record with ID %d.", job_id)
        return

    if len(pipeline) == 0:
        logger.error("Empty pipeline specified for job #%d.",
                     job_id)

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        try:
            last_result = processor(last_result)
        except Exception:
            logging.exception(("Unhandled exception caught while running processor %s in pipeline"
                               " for job #%d."),
                              processor.__name__,
                              job_id)
            last_result["success"] = False
            end_job(last_result)
        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor %s failed. Terminating pipeline for job %d.",
                         processor.__name__,
                         job_id)
            end_job(last_result)
            break

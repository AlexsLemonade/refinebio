from django.utils import timezone
from typing import List, Dict, Callable
from billiard import current_process
from data_refinery_models.models import BatchStatuses, ProcessorJob, ProcessorJobsToBatches

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
    job.worker_id = current_process().name
    job.start_time = timezone.now()
    job.save()

    batch_relations = ProcessorJobsToBatches.objects.filter(processor_job_id=job.id)
    batches = [br.batch for br in batch_relations]

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

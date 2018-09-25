from typing import List, Dict, Callable
from django.utils import timezone
from data_refinery_common.models import BatchStatuses, ProcessorJob, File
from data_refinery_common.utils import get_worker_id
from data_refinery_workers._version import __version__
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


def start_job(job_context: Dict):
    """A processor function to start jobs.

    Record in the database that this job is being started and
    retrieves the job's batches from the database and adds them to the
    dictionary passed in with the key 'batches'.
    """
    job = job_context["job"]
    job.worker_id = get_worker_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    logger.info("Starting processor Job.", processor_job=job.id)
    batches = list(job.batches.all())

    if len(batches) == 0:
        logger.error("No batches found.", processor_job=job.id)
        return {"success": False}

    job_context["batches"] = batches
    return job_context


def end_job(job_context: Dict):
    """A processor function to end jobs.

    Record in the database that this job has completed and that
    the batch has been processed if successful.
    """
    job = job_context["job"]

    if "success" in job_context:
        success = job_context["success"]
    else:
        success = True

    job.success = success
    job.end_time = timezone.now()
    job.save()

    batches = job_context["batches"] if "batches" in job_context else None

    # Clean up temp directory.
    if batches is not None and len(batches) > 0 and len(batches[0].files) > 0:
        job_dir_prefix = job_context["job_dir_prefix"] if "job_dir_prefix" in job_context else None
        batches[0].files[0].remove_temp_directory(job_dir_prefix)

    if job.success and batches is not None:
        for batch in batches:
            batch.status = BatchStatuses.PROCESSED.value
            batch.save()

        logger.info("Processor job completed successfully.", processor_job=job.id)

    # Every processor returns a dict, however end_job is always called
    # last so it doesn't need to contain anything.
    return {}


def upload_processed_files(job_context: Dict) -> Dict:
    """Uploads the processed files and removes the temp dir for the job.

    If job_context contains a "files_to_upload" key then only those
    files will be uploaded. Otherwise all files will be uploaded.
    If job_context contains a "job_dir_prefix" key then that will be
    passed through to the file methods as the `dir_name` parameter.
    """
    if "files_to_upload" in job_context:
        files = job_context["files_to_upload"]
    else:
        files = File.objects.filter(batch__in=job_context["batches"])

    if "job_dir_prefix" in job_context:
        job_dir_prefix = job_context["job_dir_prefix"]
    else:
        job_dir_prefix = None

    try:
        for file in files:
            file.upload_processed_file(job_dir_prefix)
    except Exception:
        logger.exception("Exception caught while uploading processed file %s",
                         batch=files[0].batch.id,
                         processor_job=job_context["job_id"])
        job_context["job"].failure_reason = "Exception caught while uploading processed file."
        job_context["success"] = False
        return job_context
    finally:
        # Whether or not uploading was successful, the job is over so
        # clean up the temp directory.
        files[0].remove_temp_directory(job_dir_prefix)

    return job_context


def cleanup_raw_files(job_context: Dict) -> Dict:
    """Tries to clean up raw files for the job.

    If we fail to remove the raw files, the job is still done enough
    to call a success, therefore we don't mark it as a failure.
    However logging will be important so the problem can be
    identified and the raw files cleaned up.
    """
    files = File.objects.filter(batch__in=job_context["batches"])
    for file in files:
        try:
            file.remove_raw_files()
        except:
            # If we fail to remove the raw files, the job is still done
            # enough to call a success. However logging will be important
            # so the problem can be identified and the raw files cleaned up.
            logger.exception("Exception caught while removing raw files %s",
                             file.get_temp_pre_path(),
                             batch=file.batch.id,
                             processor_job=job_context["job_id"])

    return job_context


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
        logger.error("Cannot find processor job record.", processor_job=job_id)
        return

    if len(pipeline) == 0:
        logger.error("Empty pipeline specified.",
                     procesor_job=job_id)

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        try:
            last_result = processor(last_result)
        except Exception:
            logger.exception("Unhandled exception caught while running processor %s in pipeline",
                             processor.__name__,
                             processor_job=job_id)
            last_result["success"] = False
            end_job(last_result)
        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor %s failed. Terminating pipeline.",
                         processor.__name__,
                         processor_job=job_id)
            end_job(last_result)
            break

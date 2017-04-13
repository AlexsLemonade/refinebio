from django.utils import timezone
from typing import List, Dict, Callable
from data_refinery_models.models import Batch, BatchStatuses, ProcessorJob


def start_job(kwargs: Dict):
    """Record in the database that this job is being started and
    retrieve the job and batch from the database."""
    job = (ProcessorJob
           .objects
           .filter(id=kwargs["job_id"])
           [:1]
           .get())

    job.worker_id = "For now there's only one. For now..."
    job.start_time = timezone.now()
    job.save()

    batch = (Batch
             .objects
             .filter(id=job.batch_id)
             [:1]
             .get())

    return {"job": job,
            "batch": batch}


def end_job(kwargs: Dict):
    """Record in the database that this job has completed and that
    the batch has been processed."""
    job = kwargs["job"]
    batch = kwargs["batch"]

    if "success" in kwargs:
        success = kwargs["success"]
    else:
        success = True

    job.success = success
    job.end_time = timezone.now()
    job.save()

    batch.status = BatchStatuses.PROCESSED.value
    batch.save()

    return {}


def run_pipeline(start_value: Dict, pipeline: List[Callable]):
    last_result = start_value
    for processor in pipeline:
        last_result = processor(last_result)

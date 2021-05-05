import os
import shutil

import boto3

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Pipeline, ProcessorJob, Sample
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
logger = get_and_configure_logger(__name__)

# Default to us-east-1 if the region variable can't be found
AWS_REGION = get_env_variable("AWS_REGION", "us-east-1")
AWS_BATCH_QUEUE_NAME = get_env_variable("REFINEBIO_JOB_QUEUE_NAME")

# 100 is the maximum number of jobIds you can pass to a AWS Batch's
# describe_jobs.
DESCRIBE_JOBS_PAGE_SIZE = 100

batch = boto3.client("batch", region_name=AWS_REGION)


def _find_and_remove_expired_jobs(job_context):
    """ Finds expired jobs and removes their working directories """

    job_context["deleted_items"] = []

    items_to_delete = []
    jobs_to_check = []
    for item in os.listdir(LOCAL_ROOT_DIR):

        # There may be successful processors
        if "SRP" in item or "ERP" in item or "DRP" in item:
            sub_path = os.path.join(LOCAL_ROOT_DIR, item)
            for sub_item in os.listdir(sub_path):
                try:
                    sample = Sample.objects.get(accession_code=sub_item)
                    if sample.computed_files.count() == 0:
                        # This doesn't have any associated computed files - leave it be.
                        continue
                except Sample.DoesNotExist:
                    # Interesting. This shouldn't happen at all.
                    continue
                except Exception:
                    # We can't contact the DB right now, skip deletion.
                    continue

                try:
                    sub_item_path = os.path.join(sub_path, sub_item)
                    logger.debug(
                        "Janitor deleting " + sub_item_path, contents=str(os.listdir(sub_item_path))
                    )
                    shutil.rmtree(sub_item_path)
                    job_context["deleted_items"].append(sub_item_path)
                except Exception:
                    # This job is likely vanished. No need for this directory.
                    pass

        # Processor job working directories
        if "processor_job_" in item:

            # TX Index jobs are the only ones who are allowed to hang around
            # after their jobs are finished. They're marked with an _index in their path.
            if "_index" in item:
                continue

            # Okay, does this job exist?
            try:
                job = ProcessorJob.objects.get(id=item.split("processor_job_")[1])

                jobs_to_check.append(job)
            except ProcessorJob.DoesNotExist:
                # This job has vanished from the DB - clean it up!
                logger.error("Janitor found no record of " + item + " - why?")
                items_to_delete.append(item)

    for page_start in range(0, len(jobs_to_check), DESCRIBE_JOBS_PAGE_SIZE):
        try:
            page_end = page_start + DESCRIBE_JOBS_PAGE_SIZE
            page = jobs_to_check[page_start:page_end]

            job_ids = [job.batch_job_id for job in page if job.batch_job_id]
            batch_jobs = batch.describe_jobs(jobs=job_ids)["jobs"]

            running_ids = {job["jobId"] for job in batch_jobs if job["status"] == "RUNNING"}

            for job in page:
                if not job.batch_job_id or job.batch_job_id not in running_ids:
                    # Reconstitute the path we broke apart previously.
                    items_to_delete.append("processor_job_" + str(job.id))

        except Exception:
            # We're unable to connect to the DB or Batch right
            # now, so hold onto it for right now.
            logger.exception("Problem finding job record for " + item + " - why?")
            continue

    for item in items_to_delete:
        to_delete = LOCAL_ROOT_DIR + "/" + item
        logger.debug("Janitor deleting " + to_delete)
        # This job is likely vanished. Not a problem, it's go
        shutil.rmtree(to_delete, ignore_errors=True)
        job_context["deleted_items"].append(to_delete)

    job_context["success"] = True
    return job_context


def run_janitor(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.JANITOR.value)
    job_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [utils.start_job, _find_and_remove_expired_jobs, utils.end_job],
    )
    return job_context

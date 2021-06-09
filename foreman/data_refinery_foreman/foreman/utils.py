import datetime
import sys

from django.utils import timezone

import boto3

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable

# Default to us-east-1 if the region variable can't be found
AWS_REGION = get_env_variable("AWS_REGION", "us-east-1")

logger = get_and_configure_logger(__name__)
batch = boto3.client("batch", region_name=AWS_REGION)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

PAGE_SIZE = 2000

# 100 is the maximum number of jobIds you can pass to a AWS Batch's
# describe_jobs.
DESCRIBE_JOBS_PAGE_SIZE = 100

# Setting this to a recent date will prevent the Foreman from queuing/requeuing
# jobs created before this cutoff.
JOB_CREATED_AT_CUTOFF = datetime.datetime(2019, 9, 19, tzinfo=timezone.utc)


def handle_repeated_failure(job) -> None:
    """If a job fails too many times, log it and stop retrying."""
    # Not strictly retried but will prevent the job from getting
    # retried any more times.
    job.retried = True

    # success may already be False, but if it was a hung or lost job
    # this will ensure it's marked as failed.
    job.success = False
    job.save()

    # At some point this should become more noisy/attention
    # grabbing. However for the time being just logging should be
    # sufficient because all log messages will be closely monitored
    # during early testing stages.
    logger.warn(
        "%s #%d failed %d times!!!",
        job.__class__.__name__,
        job.id,
        MAX_NUM_RETRIES + 1,
        failure_reason=job.failure_reason,
    )


def check_hung_jobs(object_list):
    hung_jobs = []
    # Batch will describe up to 100 jobs at a time.
    for page_start in range(0, len(object_list), DESCRIBE_JOBS_PAGE_SIZE):
        page_end = page_start + DESCRIBE_JOBS_PAGE_SIZE
        page = object_list[page_start:page_end]

        batch_jobs = []
        job_ids = [job.batch_job_id for job in page if job.batch_job_id]
        if job_ids:
            batch_jobs = batch.describe_jobs(jobs=job_ids)["jobs"]

        running_job_batch_ids = {job["jobId"] for job in batch_jobs if job["status"] == "RUNNING"}

        for job in page:
            if job.batch_job_id and job.batch_job_id not in running_job_batch_ids:
                hung_jobs.append(job)

    return hung_jobs


def check_lost_jobs(object_list):
    lost_jobs = []
    # Batch will describe up to 100 jobs at a time.
    for page_start in range(0, len(object_list), DESCRIBE_JOBS_PAGE_SIZE):
        page_end = page_start + DESCRIBE_JOBS_PAGE_SIZE
        page = object_list[page_start:page_end]

        batch_jobs = []
        job_ids = [job.batch_job_id for job in page if job.batch_job_id]
        if job_ids:
            batch_jobs = batch.describe_jobs(jobs=job_ids)["jobs"]

        # Need to ignore statuses where the job wouldn't have its
        # start_time set. This includes RUNNING because it may not
        # have yet gotten to that point.
        ignore = ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]
        ignore_job_batch_ids = {job["jobId"] for job in batch_jobs if job["status"] in ignore}

        for job in page:
            if not job.batch_job_id or job.batch_job_id not in ignore_job_batch_ids:
                lost_jobs.append(job)

    return lost_jobs

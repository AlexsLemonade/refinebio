import datetime
import sys

from django.utils import timezone

import boto3

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable, get_env_variable_gracefully

# Default to us-east-1 if the region variable can't be found
AWS_REGION = get_env_variable("AWS_REGION", "us-east-1")
AWS_BATCH_QUEUE_NAME = get_env_variable("REFINEBIO_JOB_QUEUE_NAME")

logger = get_and_configure_logger(__name__)
batch = boto3.client("batch", region_name=AWS_REGION)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

# This can be overritten by the env var "MAX_TOTAL_JOBS"
DEFAULT_MAX_JOBS = 5000

# This is the absolute max number of downloader jobs that should ever
# be queued across the whole cluster no matter how many nodes we
# have. This is important because too many downloader jobs and we take
# down NCBI.
HARD_MAX_DOWNLOADER_JOBS = 1000

PAGE_SIZE = 2000

# 100 is the maximum number of jobIds you can pass to a AWS Batch's
# describe_jobs.
DESCRIBE_JOBS_PAGE_SIZE = 100

# Setting this to a recent date will prevent the Foreman from queuing/requeuing
# jobs created before this cutoff.
JOB_CREATED_AT_CUTOFF = datetime.datetime(2019, 9, 19, tzinfo=timezone.utc)

# The total number of jobs currently in the Batch queue
JOBS_IN_QUEUE = 0
# The number of downloader jobs currently in the Batch queue
DOWNLOADER_JOBS_IN_QUEUE = 0
TIME_OF_LAST_JOB_CHECK = timezone.now() - datetime.timedelta(minutes=10)


def increment_job_queue():
    global JOBS_IN_QUEUE
    JOBS_IN_QUEUE += 1


def increment_downloader_job_queue():
    global DOWNLOADER_JOBS_IN_QUEUE
    DOWNLOADER_JOBS_IN_QUEUE += 1


def count_jobs_in_queue(window=datetime.timedelta(minutes=2)) -> int:
    """Counts how many jobs are in the job queue that aren't finished.

    Also returns how many of those are downloader jobs."""
    global JOBS_IN_QUEUE
    global DOWNLOADER_JOBS_IN_QUEUE
    global TIME_OF_LAST_JOB_CHECK

    if timezone.now() - TIME_OF_LAST_JOB_CHECK > window:
        try:
            num_jobs = 0
            num_downloader_jobs = 0

            # AWS Batch only returns one status at a time and doesn't provide a `count` or `total`.
            for status in ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]:
                list_jobs_dict = batch.list_jobs(jobQueue=AWS_BATCH_QUEUE_NAME, jobStatus=status)

                def is_downloader(job_json):
                    return job_json["jobName"].startswith("Downloader")

                num_jobs += len(list_jobs_dict["jobSummaryList"])
                num_downloader_jobs += len(
                    list(filter(is_downloader, list_jobs_dict["jobSummaryList"]))
                )
                while "nextToken" in list_jobs_dict and list_jobs_dict["nextToken"]:
                    list_jobs_dict = batch.list_jobs(
                        jobQueue=AWS_BATCH_QUEUE_NAME,
                        jobStatus=status,
                        nextToken=list_jobs_dict["nextToken"],
                    )
                    num_jobs += len(list_jobs_dict["jobSummaryList"])
                    num_downloader_jobs += len(
                        list(filter(is_downloader, list_jobs_dict["jobSummaryList"]))
                    )

            JOBS_IN_QUEUE = num_jobs
            DOWNLOADER_JOBS_IN_QUEUE = num_downloader_jobs
        except Exception:
            logger.exception("Unable to determine number of Batch jobs.")
            # Can't query Batch, return an impossibly high number to prevent
            # additonal queuing from happening:
            JOBS_IN_QUEUE = sys.maxsize
            DOWNLOADER_JOBS_IN_QUEUE = sys.maxsize

        TIME_OF_LAST_JOB_CHECK = timezone.now()

    return JOBS_IN_QUEUE, DOWNLOADER_JOBS_IN_QUEUE


def get_capacity_for_jobs() -> bool:
    """Returns how many jobs the queue has capacity for.
    """
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = int(get_env_variable_gracefully("MAX_TOTAL_JOBS", DEFAULT_MAX_JOBS))

    num_current_jobs, _ = count_jobs_in_queue()

    return max(MAX_TOTAL_JOBS - num_current_jobs, 0)


def get_capacity_for_downloader_jobs() -> int:
    """Returns how many downloader jobs the queue has capacity for.

    This has to respect the overall job limit and also the downloader job limit.
    """
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = int(get_env_variable_gracefully("MAX_TOTAL_JOBS", DEFAULT_MAX_JOBS))

    num_jobs, num_downloader_jobs = count_jobs_in_queue()

    downloader_capacity = HARD_MAX_DOWNLOADER_JOBS - num_downloader_jobs
    overall_capacity = MAX_TOTAL_JOBS - num_jobs

    if downloader_capacity < 0 or overall_capacity < 0:
        return 0

    return min(downloader_capacity, overall_capacity)


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

        job_ids = [job.batch_job_id for job in page if job.batch_job_id]
        batch_jobs = batch.describe_jobs(jobs=job_ids)["jobSummaryList"]

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

        job_ids = [job.batch_job_id for job in page if job.batch_job_id]
        batch_jobs = batch.describe_jobs(jobs=job_ids)["jobSummaryList"]

        # Need to ignore statuses where the job wouldn't have its
        # start_time set. This includes RUNNING because it may not
        # have yet gotten to that point.
        ignore = ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]
        ignore_job_batch_ids = {job["jobId"] for job in batch_jobs if job["status"] in ignore}

        for job in page:
            if not job.batch_job_id or job.batch_job_id not in ignore_job_batch_ids:
                lost_jobs.append(job)

    return lost_jobs

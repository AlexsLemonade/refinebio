"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals

import datetime
import sys
from enum import Enum

from django.conf import settings
from django.utils import timezone

import boto3

from data_refinery_common.job_lookup import (
    SMASHER_JOB_TYPES,
    Downloaders,
    ProcessorPipeline,
    SurveyJobTypes,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)


AWS_REGION = get_env_variable(
    "AWS_REGION", "us-east-1"
)  # Default to us-east-1 if the region variable can't be found

# Job definitons are AWS objects so they have to be namespaced for our stack.
JOB_DEFINITION_PREFIX = get_env_variable("JOB_DEFINITION_PREFIX", "")

# These two constants refer to image names that should be used for
# multiple jobs.
BATCH_TRANSCRIPTOME_JOB = "TRANSCRIPTOME_INDEX"
BATCH_DOWNLOADER_JOB = "DOWNLOADER"
NONE_JOB_ERROR_TEMPLATE = "send_job was called with NONE job_type: {} for {} job {}"
TIME_OF_LAST_JOB_CHECK = timezone.now() - datetime.timedelta(minutes=10)


batch = boto3.client("batch", region_name=AWS_REGION)

# Initialize the queue depths to zero to be used as a global.
JOB_QUEUE_DEPTHS = {}
for queue_name in settings.AWS_BATCH_QUEUE_ALL_NAMES:
    JOB_QUEUE_DEPTHS[queue_name] = 0

DOWNLOADER_JOB_QUEUE_DEPTHS = {}
for queue_name in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
    DOWNLOADER_JOB_QUEUE_DEPTHS[queue_name] = 0


def count_jobs_in_queue(batch_job_queue) -> dict:
    """Counts how many jobs are in the job queue that aren't finished.

    Also returns how many of those are downloader jobs."""
    try:
        num_jobs = 0
        num_downloader_jobs = 0

        # AWS Batch only returns one status at a time and doesn't provide a `count` or `total`.
        for status in ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]:
            list_jobs_dict = batch.list_jobs(jobQueue=batch_job_queue, jobStatus=status)

            def is_downloader(job_json):
                return job_json["jobName"].startswith("Downloader")

            num_jobs += len(list_jobs_dict["jobSummaryList"])
            num_downloader_jobs += len(
                list(filter(is_downloader, list_jobs_dict["jobSummaryList"]))
            )
            while "nextToken" in list_jobs_dict and list_jobs_dict["nextToken"]:
                list_jobs_dict = batch.list_jobs(
                    jobQueue=batch_job_queue,
                    jobStatus=status,
                    nextToken=list_jobs_dict["nextToken"],
                )
                num_jobs += len(list_jobs_dict["jobSummaryList"])
                num_downloader_jobs += len(
                    list(filter(is_downloader, list_jobs_dict["jobSummaryList"]))
                )
    except Exception:
        logger.exception("Unable to determine number of Batch jobs.")
        # Can't query Batch, return an impossibly high number to prevent
        # additonal queuing from happening:
        num_jobs = sys.maxsize
        num_downloader_jobs = sys.maxsize

    return {"all_jobs": num_jobs, "downloader_jobs": num_downloader_jobs}


def get_job_queue_depths(window=datetime.timedelta(minutes=2)):
    global TIME_OF_LAST_JOB_CHECK
    global JOB_QUEUE_DEPTHS
    global DOWNLOADER_JOB_QUEUE_DEPTHS

    if timezone.now() - TIME_OF_LAST_JOB_CHECK > window:
        for job_queue in settings.AWS_BATCH_QUEUE_ALL_NAMES:
            job_count_dict = count_jobs_in_queue(job_queue)
            JOB_QUEUE_DEPTHS[job_queue] = job_count_dict["all_jobs"]

            if job_queue in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
                DOWNLOADER_JOB_QUEUE_DEPTHS[job_queue] = job_count_dict["downloader_jobs"]

        TIME_OF_LAST_JOB_CHECK = timezone.now()

    return {"all_jobs": JOB_QUEUE_DEPTHS, "downloader_jobs": DOWNLOADER_JOB_QUEUE_DEPTHS}


def get_job_queue_depth(job_queue_name):
    return get_job_queue_depths()["all_jobs"][job_queue_name]


def get_downloader_job_queue_depth(job_queue_name):
    return get_job_queue_depths()["downloader_jobs"][job_queue_name]


def get_capacity_for_jobs() -> bool:
    """Returns how many jobs the system has remaining capacity for.
    """
    num_job_queues = len(settings.AWS_BATCH_QUEUE_WORKERS_NAMES)
    MAX_JOBS_PER_NODE = int(get_env_variable("MAX_JOBS_PER_NODE"))
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = MAX_JOBS_PER_NODE * num_job_queues

    job_queue_depths = get_job_queue_depths()["all_jobs"]
    num_current_jobs = 0
    for job_queue in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
        num_current_jobs += job_queue_depths[job_queue]

    return max(MAX_TOTAL_JOBS - num_current_jobs, 0)


def get_capacity_for_downloader_jobs() -> int:
    """Returns how many downloader jobs the system has remaining capacity for.

    This has to respect the overall job limit and also the downloader job limit.
    """
    num_job_queues = len(settings.AWS_BATCH_QUEUE_WORKERS_NAMES)
    MAX_DOWNLOADER_JOBS_PER_NODE = settings.MAX_DOWNLOADER_JOBS_PER_NODE
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_DOWNLOADER_JOBS = MAX_DOWNLOADER_JOBS_PER_NODE * num_job_queues

    job_queue_depths = get_job_queue_depths()["downloader_jobs"]
    num_downloader_jobs = 0
    for job_queue in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
        num_downloader_jobs += job_queue_depths[job_queue]

    downloader_capacity = MAX_TOTAL_DOWNLOADER_JOBS - num_downloader_jobs
    overall_capacity = get_capacity_for_jobs()

    if downloader_capacity < 0 or overall_capacity < 0:
        return 0

    return min(downloader_capacity, overall_capacity)


def increment_job_queue_depth(job_queue_name):
    global DOWNLOADER_JOB_QUEUE_DEPTHS
    JOB_QUEUE_DEPTHS[job_queue_name] += 1


def increment_downloader_job_queue_depth(job_queue_name):
    global DOWNLOADER_JOB_QUEUE_DEPTHS
    DOWNLOADER_JOB_QUEUE_DEPTHS[job_queue_name] += 1


def get_batch_queue_for_downloader_job():
    """Logic for distributing downloader jobs across queues.
    """
    job_queue_depths = get_job_queue_depths()["downloader_jobs"]
    for job_queue in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
        if job_queue_depths[job_queue] <= settings.MAX_DOWNLOADER_JOBS_PER_NODE:
            return job_queue

    # If none of the job queues have capacity for downloader
    # jobs do not queue the job. This ensures we actually
    # process data we download rather than overloading the
    # instance with downloader jobs until its disk is full.
    return None


def get_first_job_queue_with_capacity():
    """Returns the first job queue that has capacity for more jobs.

    If there are no job queues with capacity, returns None.
    """
    job_queue_depths = get_job_queue_depths()["all_jobs"]
    for job_queue in settings.AWS_BATCH_QUEUE_WORKERS_NAMES:
        if job_queue_depths[job_queue] <= settings.MAX_JOBS_PER_NODE:
            return job_queue

    return None


def get_batch_queue_for_job(job_type, job):
    if job_type is ProcessorPipeline.SMASHER:
        return settings.AWS_BATCH_QUEUE_SMASHER_NAME
    elif job_type in [ProcessorPipeline.CREATE_COMPENDIA, ProcessorPipeline.CREATE_QUANTPENDIA]:
        # This should probably use the workers' queue for everything
        # but human/mouse.
        # https://github.com/AlexsLemonade/refinebio/issues/2744
        return settings.AWS_BATCH_QUEUE_COMPENDIA_NAME
    elif job_type in list(Downloaders):
        return get_batch_queue_for_downloader_job()
    elif job_type in list(ProcessorPipeline):
        log_str = f"The downloader job is {job.downloader_job}."
        if job.downloader_job:
            log_str = log_str + f", it's queue is {job.downloader_job.batch_job_queue}"

        logger.info(log_str, job_type=job_type)
        # Queue it in the same queue as the downloader job as long as
        # that is set and still available.
        if (
            not job.downloader_job
            or not job.downloader_job.batch_job_queue
            or job.downloader_job.batch_job_queue not in settings.AWS_BATCH_QUEUE_WORKERS_NAMES
        ):
            return get_first_job_queue_with_capacity()
        else:
            return job.downloader_job.batch_job_queue

    elif job_type in list(SurveyJobTypes):
        # We always want to queue a survey job, so just look for the queue
        # with the smallest number of jobs in it. The only time we return
        # None is if there's no job queues for some reason.
        return get_first_job_queue_with_capacity()
    else:
        # Handle the case where it's none of the above. Shouldn't happen.
        raise ValueError(f"Job id {job.id} had an invalid job_type: {job_type.value}")


def get_job_name(job_type, job_id):
    if (
        job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG
        or job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT
    ):
        return BATCH_TRANSCRIPTOME_JOB
    elif job_type is ProcessorPipeline.SALMON:
        return ProcessorPipeline.SALMON.value
    elif job_type is ProcessorPipeline.TXIMPORT:
        return ProcessorPipeline.TXIMPORT.value
    elif job_type is ProcessorPipeline.AFFY_TO_PCL:
        return ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type is ProcessorPipeline.NO_OP:
        return ProcessorPipeline.NO_OP.value
    elif job_type is ProcessorPipeline.ILLUMINA_TO_PCL:
        return ProcessorPipeline.ILLUMINA_TO_PCL.value
    elif job_type is ProcessorPipeline.SMASHER:
        return ProcessorPipeline.SMASHER.value
    elif job_type is ProcessorPipeline.JANITOR:
        return ProcessorPipeline.JANITOR.value
    elif job_type is ProcessorPipeline.QN_REFERENCE:
        return ProcessorPipeline.QN_REFERENCE.value
    elif job_type is ProcessorPipeline.CREATE_COMPENDIA:
        return ProcessorPipeline.CREATE_COMPENDIA.value
    elif job_type is ProcessorPipeline.CREATE_QUANTPENDIA:
        return ProcessorPipeline.CREATE_QUANTPENDIA.value
    elif job_type is ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL:
        # Agilent twocolor uses the same job specification as Affy.
        return ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type in list(Downloaders):
        return BATCH_DOWNLOADER_JOB
    elif job_type in list(SurveyJobTypes):
        return job_type.value
    elif job_type is Downloaders.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Downloader", job_id))
    elif job_type is ProcessorPipeline.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Processor", job_id))
    else:
        raise ValueError("Invalid job_type: {}".format(job_type.value))


def is_job_processor(job_type):
    return job_type not in list(Downloaders) and job_type not in list(SurveyJobTypes)


def send_job(job_type: Enum, job, is_dispatch=False) -> bool:
    job_name = get_job_name(job_type, job.id)
    is_processor = is_job_processor(job_type)

    if settings.AUTO_DISPATCH_BATCH_JOBS:
        # We only want to dispatch processor jobs directly.
        # Everything else will be handled by the Foreman, which will increment the retry counter.
        should_dispatch = is_processor or is_dispatch or (not settings.RUNNING_IN_CLOUD)
    else:
        should_dispatch = is_dispatch  # only dispatch when specifically requested to

    if should_dispatch:
        batch = boto3.client("batch", region_name=AWS_REGION)

        job_name = JOB_DEFINITION_PREFIX + job_name

        # Smasher related and tximport jobs  don't have RAM tiers.
        if job_type not in SMASHER_JOB_TYPES and job_type is not ProcessorPipeline.TXIMPORT:
            job_name = job_name + "_" + str(job.ram_amount)

        job_queue = get_batch_queue_for_job(job_type, job)

        try:
            batch_response = batch.submit_job(
                jobName=job_name + f"_{job.id}",
                jobQueue=job_queue,
                jobDefinition=job_name,
                parameters={"job_name": job_type.value, "job_id": str(job.id)},
            )
            job.batch_job_queue = job_queue
            job.batch_job_id = batch_response["jobId"]
            job.save()

            increment_job_queue_depth(job_queue)
            if job_type in list(Downloaders):
                increment_downloader_job_queue_depth(job_queue)

            return True
        except Exception as e:
            logger.warn(
                "Unable to Dispatch Batch Job.",
                job_name=job_type.value,
                job_id=str(job.id),
                reason=str(e),
            )
            raise
    else:
        job.num_retries = job.num_retries - 1
        job.save()

    return True

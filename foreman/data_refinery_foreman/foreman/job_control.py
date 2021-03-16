import datetime
import sys
import time
from typing import List

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
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    ComputedFile,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
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
# This number is related to the size of the nomad lead server instance.
# The larger the instance, the more jobs it can store in its memory at once.
# We've found these limits by testing:
#  * t2.medium can handle 5000 jobs.
#  * m5.xlarge can hanlde 20000 jobs.
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

# The number of downloader jobs currently in the nomad queue
JOBS_IN_QUEUE = 0
DOWNLOADER_JOBS_IN_QUEUE = 0
TIME_OF_LAST_JOB_CHECK = timezone.now() - datetime.timedelta(minutes=10)

# The minimum amount of time in between each iteration of the main
# loop. We could loop much less frequently than every two minutes if
# the work we do takes longer than 2 minutes, but this will prevent
# excessive spinning.
MIN_LOOP_TIME = datetime.timedelta(seconds=15)

# How frequently we dispatch Janitor jobs and clean unplaceable jobs
# out of the Nomad queue.
JANITOR_DISPATCH_TIME = datetime.timedelta(minutes=30)

# How frequently we clean up the database.
DBCLEAN_TIME = datetime.timedelta(hours=6)

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


def requeue_downloader_job(last_job: DownloaderJob) -> (bool, str):
    """Queues a new downloader job.

    The new downloader job will have num_retries one greater than
    last_job.num_retries.

    Returns True and the volume index of the downloader job upon successful dispatching,
    False and an empty string otherwise.
    """
    num_retries = last_job.num_retries + 1

    ram_amount = last_job.ram_amount
    # If there's no start time then it's likely that the instance got
    # cycled which means we didn't get OOM-killed, so we don't need to
    # increase the RAM amount.
    if last_job.start_time and last_job.failure_reason is None:
        if ram_amount == 1024:
            ram_amount = 4096
        elif ram_amount == 4096:
            ram_amount = 16384

    original_file = last_job.original_files.first()

    if not original_file:
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = (
            "Foreman told to requeue a DownloaderJob without an OriginalFile - why?!"
        )
        last_job.save()
        logger.info(
            "Foreman told to requeue a DownloaderJob without an OriginalFile - why?!",
            last_job=str(last_job),
        )
        return False, ""

    if not original_file.needs_processing():
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = "Foreman told to redownload job with prior successful processing."
        last_job.save()
        logger.info(
            "Foreman told to redownload job with prior successful processing.",
            last_job=str(last_job),
        )
        return False, ""

    first_sample = original_file.samples.first()

    # This is a magic string that all the dbGaP studies appear to have
    if first_sample and ("in the dbGaP study" in first_sample.title):
        last_job.no_retry = True
        last_job.success = False
        last_job.failure_reason = "Sample is dbGaP access controlled."
        last_job.save()
        logger.info(
            "Avoiding requeuing for DownloaderJob for dbGaP run accession: "
            + str(first_sample.accession_code)
        )
        return False, ""

    new_job = DownloaderJob(
        num_retries=num_retries,
        downloader_task=last_job.downloader_task,
        ram_amount=ram_amount,
        accession_code=last_job.accession_code,
        was_recreated=last_job.was_recreated,
        # For now default to 0, will need to be calculated.
        volume_index=0,
    )
    new_job.save()

    for original_file in last_job.original_files.all():
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=new_job, original_file=original_file
        )

    logger.debug(
        "Requeuing Downloader Job which had ID %d with a new Downloader Job with ID %d.",
        last_job.id,
        new_job.id,
    )
    try:
        if send_job(Downloaders[last_job.downloader_task], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
            return False, ""
    except Exception:
        logger.error(
            "Failed to requeue DownloaderJob which had ID %d with a new DownloaderJob with ID %d.",
            last_job.id,
            new_job.id,
        )
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()
        return False, ""

    return True, new_job.volume_index


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


def handle_downloader_jobs(jobs: List[DownloaderJob]) -> None:
    """For each job in jobs, either retry it or log it.

    No more than queue_capacity jobs will be retried.
    """
    global DOWNLOADER_JOBS_IN_QUEUE

    queue_capacity = get_capacity_for_downloader_jobs()

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if jobs_dispatched >= queue_capacity:
            logger.info(
                "We hit the maximum downloader jobs / capacity ceiling, "
                "so we're not handling any more downloader jobs now."
            )
            return

        if job.num_retries < MAX_NUM_RETRIES:
            requeue_success, dispatched_volume = requeue_downloader_job(job)
            if requeue_success:
                jobs_dispatched = jobs_dispatched + 1
                DOWNLOADER_JOBS_IN_QUEUE += 1
        else:
            handle_repeated_failure(job)


def retry_failed_downloader_jobs() -> None:
    """Handle downloader jobs that were marked as a failure."""
    failed_jobs = (
        DownloaderJob.failed_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF)
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )

    paginator = Paginator(failed_jobs, PAGE_SIZE, "created_at")
    page = paginator.page()
    page_count = 0

    if len(page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_downloader_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling failed (explicitly-marked-as-failure) downloader jobs "
            "because there is no capacity for them."
        )

    while queue_capacity > 0:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) downloader jobs!", page_count
        )

        handle_downloader_jobs(page.object_list)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
            queue_capacity = get_capacity_for_downloader_jobs()
        else:
            break


def check_hung_jobs(object_list):
    hung_jobs = []
    # Batch will describe up to 100 jobs at a time.
    for page_start in range(0, len(object_list), DESCRIBE_JOBS_PAGE_SIZE):
        page_end = page_start + DESCRIBE_JOBS_PAGE_SIZE
        network_page = object_list[page_start:page_end]

        job_ids = [job.batch_job_id for job in network_page if job.batch_job_id]
        batch_jobs = batch.describe_jobs(jobs=job_ids)["jobs"]

        hung_job_batch_ids = {job["jobId"] for job in batch_jobs if job["status"] != "RUNNING"}

        for job in network_page:
            if job.batch_job_id and job.batch_job_id in hung_job_batch_ids:
                hung_jobs.append(job)

    return hung_jobs


def retry_hung_downloader_jobs() -> None:
    """Retry downloader jobs that were started but never finished."""
    potentially_hung_jobs = (
        DownloaderJob.hung_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF)
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )
    paginator = Paginator(potentially_hung_jobs, PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_downloader_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling failed (explicitly-marked-as-failure) downloader jobs "
            "because there is no capacity for them."
        )
    while queue_capacity > 0:
        hung_jobs = check_hung_jobs(database_page.object_list)

        if hung_jobs:
            logger.info(
                "Handling page %d of hung (started-but-never-finished) downloader jobs!",
                database_page_count,
                jobs_count=len(hung_jobs),
            )
            handle_downloader_jobs(hung_jobs)

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = get_capacity_for_downloader_jobs()
        else:
            break


def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    # The Salmon pipeline is quite RAM-sensitive.
    # Try it again with an increased RAM amount, if possible.
    new_ram_amount = last_job.ram_amount

    # If there's no start time then it's likely that the instance got
    # cycled which means we didn't get OOM-killed, so we don't need to
    # increase the RAM amount.
    if last_job.start_time:
        # These initial values are set in common/job_lookup.py:determine_ram_amount
        if (
            last_job.pipeline_applied == "SALMON"
            or last_job.pipeline_applied == "TXIMPORT"
            or last_job.pipeline_applied.startswith("TRANSCRIPTOME")
        ):
            if new_ram_amount == 4096:
                new_ram_amount = 8192
            if new_ram_amount == 8192:
                new_ram_amount = 12288
            elif new_ram_amount == 12288:
                new_ram_amount = 16384
            elif new_ram_amount == 16384:
                new_ram_amount = 32768
            elif new_ram_amount == 32768:
                new_ram_amount = 65536
        # The AFFY pipeline is somewhat RAM-sensitive.
        # Also NO_OP can fail and be retried, so we want to attempt ramping up ram.
        # Try it again with an increased RAM amount, if possible.
        elif last_job.pipeline_applied == "AFFY_TO_PCL" or last_job.pipeline_applied == "NO_OP":
            if new_ram_amount == 2048:
                new_ram_amount = 4096
            elif new_ram_amount == 4096:
                new_ram_amount = 8192
            elif new_ram_amount == 8192:
                new_ram_amount = 32768
        elif (
            last_job.pipeline_applied == "ILLUMINA_TO_PCL"
            and "non-zero exit status -9" in last_job.failure_reason
        ):
            if new_ram_amount == 2048:
                new_ram_amount = 4096
            elif new_ram_amount == 4096:
                new_ram_amount = 8192

    volume_index = last_job.volume_index
    # Make sure volume_index is set to something, unless it's a
    # smasher job type because the smasher instance doesn't have a
    # volume_index.
    if (not volume_index or volume_index == "-1") and ProcessorPipeline[
        last_job.pipeline_applied
    ] not in SMASHER_JOB_TYPES:
        # Default to 0 for now. At some point we'll have to detect
        # when 0 is full, and if so go to 1, etc.
        volume_index = 0

    new_job = ProcessorJob(
        num_retries=num_retries,
        pipeline_applied=last_job.pipeline_applied,
        ram_amount=new_ram_amount,
        volume_index=volume_index,
    )
    new_job.save()

    for original_file in last_job.original_files.all():
        ProcessorJobOriginalFileAssociation.objects.get_or_create(
            processor_job=new_job, original_file=original_file
        )

    for dataset in last_job.datasets.all():
        ProcessorJobDatasetAssociation.objects.get_or_create(processor_job=new_job, dataset=dataset)

    try:
        logger.debug(
            "Requeuing Processor Job which had ID %d with a new Processor Job with ID %d.",
            last_job.id,
            new_job.id,
        )
        if send_job(ProcessorPipeline[last_job.pipeline_applied], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
    except Exception:
        logger.warn(
            "Failed to requeue Processor Job which had ID %d with a new Processor Job with ID %d.",
            last_job.id,
            new_job.id,
            exc_info=1,
        )
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()


def handle_processor_jobs(
    jobs: List[ProcessorJob], queue_capacity: int = None, ignore_ceiling=False
) -> None:
    """For each job in jobs, either retry it or log it.

    No more than queue_capacity jobs will be retried.
    """
    if queue_capacity is None:
        queue_capacity = get_capacity_for_jobs()

    jobs_dispatched = 0
    for count, job in enumerate(jobs):

        if not ignore_ceiling and jobs_dispatched >= queue_capacity:
            logger.info(
                "We hit the maximum total jobs ceiling, "
                "so we're not handling any more processor jobs now."
            )
            return

        if job.num_retries < MAX_NUM_RETRIES:
            requeue_processor_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)


def retry_failed_processor_jobs() -> None:
    """Handle processor jobs that were marked as a failure.

    Ignores Janitor jobs since they are queued every half hour anyway."""
    failed_jobs = (
        ProcessorJob.failed_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF)
        .exclude(pipeline_applied="JANITOR")
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )

    paginator = Paginator(failed_jobs, 200, "created_at")
    page = paginator.page()
    page_count = 0

    if len(page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    if queue_capacity > 0:
        for i in range(queue_capacity):
            logger.info(
                "Handling page %d of failed (explicitly-marked-as-failure) processor jobs!",
                page_count,
            )
            handle_processor_jobs(page.object_list, queue_capacity)

            if page.has_next():
                page = paginator.page(page.next_page_number())
                page_count = page_count + 1
                queue_capacity = get_capacity_for_jobs()
            else:
                break


def retry_hung_processor_jobs() -> None:
    """Retry processor jobs that were started but never finished."""
    potentially_hung_jobs = (
        ProcessorJob.hung_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF)
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )
    paginator = Paginator(potentially_hung_jobs, PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling failed (explicitly-marked-as-failure) processor jobs "
            "because there is no capacity for them."
        )
    while queue_capacity > 0:
        hung_jobs = check_hung_jobs(database_page.object_list)

        if hung_jobs:
            logger.info(
                "Handling page %d of hung (started-but-never-finished) processor jobs!",
                database_page_count,
                jobs_count=len(hung_jobs),
            )
            handle_processor_jobs(hung_jobs)

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = get_capacity_for_jobs()
        else:
            break


def requeue_survey_job(last_job: SurveyJob) -> None:
    """Queues a new survey job.

    The new survey job will have num_retries one greater than
    last_job.num_retries.
    """

    num_retries = last_job.num_retries + 1

    new_job = SurveyJob(num_retries=num_retries, source_type=last_job.source_type)

    if new_job.num_retries == 1:
        new_job.ram_amount = 4096
    elif new_job.num_retries in [2, 3]:
        new_job.ram_amount = 16384
    else:
        new_job.ram_amount = 256

    new_job.save()

    keyvalues = SurveyJobKeyValue.objects.filter(survey_job=last_job)

    for keyvalue in keyvalues:
        SurveyJobKeyValue.objects.get_or_create(
            survey_job=new_job, key=keyvalue.key, value=keyvalue.value,
        )

    logger.debug(
        "Requeuing SurveyJob which had ID %d with a new SurveyJob with ID %d.",
        last_job.id,
        new_job.id,
    )

    try:
        if send_job(SurveyJobTypes.SURVEYOR, job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            new_job.delete()
    except Exception:
        logger.error(
            "Failed to requeue Survey Job which had ID %d with a new Surevey Job with ID %d.",
            last_job.id,
            new_job.id,
        )
        # Can't communicate with AWS just now, leave the job for a later loop.
        new_job.delete()

    return True


def handle_survey_jobs(jobs: List[SurveyJob], queue_capacity: int = None) -> None:
    """For each job in jobs, either retry it or log it.

    No more than queue_capacity jobs will be retried.
    """
    if queue_capacity is None:
        queue_capacity = get_capacity_for_jobs()

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if jobs_dispatched >= queue_capacity:
            logger.info(
                "We hit the maximum total jobs ceiling,"
                " so we're not handling any more survey jobs now."
            )
            return

        if job.num_retries < MAX_NUM_RETRIES:
            requeue_survey_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)


def retry_failed_survey_jobs() -> None:
    """Handle survey jobs that were marked as a failure."""
    failed_jobs = SurveyJob.failed_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).order_by(
        "pk"
    )

    paginator = Paginator(failed_jobs, 200)
    page = paginator.page()
    page_count = 0

    if len(page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    while queue_capacity > 0:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) survey jobs!", page_count
        )
        handle_survey_jobs(page.object_list, queue_capacity)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
            queue_capacity = get_capacity_for_jobs()
        else:
            break


##
# Janitor
##


def send_janitor_jobs():
    """Dispatch a Janitor job for each compute environment.

    For now that's just one, but this will eventually need to loop.
    """
    new_job = ProcessorJob(num_retries=0, pipeline_applied="JANITOR", ram_amount=2048)
    new_job.save()
    logger.info("Sending Janitor Job.", job_id=new_job.id)
    try:
        send_job(ProcessorPipeline["JANITOR"], job=new_job, is_dispatch=True)
    except Exception:
        # If we can't dispatch this job, something else has gone wrong, we can get it next loop.
        return


def clean_database():
    """ Removes duplicated objects that may have appeared through race, OOM, bugs, etc.
    See: https://github.com/AlexsLemonade/refinebio/issues/1183
    """

    # Hide smashable files
    computed_files = ComputedFile.objects.filter(s3_bucket=None, s3_key=None, is_smashable=True)
    logger.info("Cleaning unsynced files!", num_to_clean=computed_files.count())

    # We don't do this in bulk because we want the properties set by save() as well
    for computed_file in computed_files:
        computed_file.is_public = False
        computed_file.save()

    logger.info("Cleaned files!")


##
# Main loop
##


def monitor_jobs():
    """Main Foreman thread that helps manage the Nomad job queue.

    Will find jobs that failed, hung, or got lost and requeue them.

    Also will queue up Janitor jobs regularly to free up disk space.

    Also cleans jobs out of the Nomad queue which cannot be queued
    because the volume containing the job's data isn't mounted.

    It does so on a loop forever that won't spin faster than
    MIN_LOOP_TIME, but it may spin slower than that.
    """
    last_janitorial_time = timezone.now()
    last_dbclean_time = timezone.now()

    while True:
        # Perform two heartbeats, one for the logs and one for Monit:
        logger.info("The Foreman's heart is beating, but he does not feel.")

        # Write the health file for Monit to check
        now_secs = int(time.time())
        with open("/tmp/foreman_last_time", "w") as timefile:
            timefile.write(str(now_secs))

        start_time = timezone.now()

        # Requeue jobs of each failure class for each job type.
        # The order of processor -> downloader -> surveyor is intentional.
        # Processors go first so we process data sitting on disk.
        # Downloaders go first so we actually queue up the jobs in the database.
        # Surveyors go last so we don't end up with tons and tons of unqueued jobs.
        requeuing_functions_in_order = [
            retry_failed_processor_jobs,
            # retry_hung_processor_jobs,
            # retry_lost_processor_jobs,
            retry_failed_downloader_jobs,
            # retry_hung_downloader_jobs,
            # retry_lost_downloader_jobs,
            retry_failed_survey_jobs,
            # retry_hung_survey_jobs,
            # retry_lost_survey_jobs,
        ]

        for function in requeuing_functions_in_order:
            try:
                function()
            except Exception:
                logger.exception("Caught exception in %s: ", function.__name__)

        if settings.RUNNING_IN_CLOUD:
            if timezone.now() - last_janitorial_time > JANITOR_DISPATCH_TIME:
                send_janitor_jobs()
                last_janitorial_time = timezone.now()

            if timezone.now() - last_dbclean_time > DBCLEAN_TIME:
                clean_database()
                last_dbclean_time = timezone.now()

        loop_time = timezone.now() - start_time
        if loop_time < MIN_LOOP_TIME:
            remaining_time = MIN_LOOP_TIME - loop_time
            if remaining_time.seconds > 0:
                time.sleep(remaining_time.seconds)

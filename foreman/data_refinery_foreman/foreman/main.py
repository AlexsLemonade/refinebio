import time
import nomad
from nomad.api.exceptions import URLNotFoundNomadException
from typing import Callable, List
from threading import Thread
from functools import wraps
from retrying import retry
from datetime import timedelta
from django.utils import timezone
from django.db import transaction
from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation,
    ProcessorJobDatasetAssociation
)
from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable


logger = get_and_configure_logger(__name__)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

# For now it seems like no jobs should take longer than a day to be
# either picked up or run.
MAX_DOWNLOADER_RUN_TIME = timedelta(days=1)
MAX_PROCESSOR_RUN_TIME = timedelta(days=1)
MAX_QUEUE_TIME = timedelta(days=1)

# To prevent excessive spinning loop no more than once every 30 minutes.
MIN_LOOP_TIME = timedelta(minutes=30)
THREAD_WAIT_TIME = 10.0


@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_downloader_job(last_job: DownloaderJob) -> None:
    """Queues a new downloader job.

    The new downloader job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    new_job = DownloaderJob(num_retries=num_retries,
                            downloader_task=last_job.downloader_task,
                            accession_code=last_job.accession_code)
    new_job.save()

    for original_file in last_job.original_files.all():
        DownloaderJobOriginalFileAssociation.objects.get_or_create(downloader_job=new_job,
                                                           original_file=original_file)

    logger.info("Requeuing Downloader Job which had ID %d with a new Downloader Job with ID %d.",
                last_job.id,
                new_job.id)
    send_job(Downloaders[last_job.downloader_task], new_job)

    last_job.retried = True
    last_job.success = False
    last_job.retried_job = new_job
    last_job.save()


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
    logger.warn("%s #%d failed %d times!!!", job.__class__.__name__, job.id, MAX_NUM_RETRIES + 1)


def handle_downloader_jobs(jobs: List[DownloaderJob]) -> None:
    """For each job in jobs, either retry it or log it."""
    for job in jobs:
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_downloader_job(job)
        else:
            handle_repeated_failure(job)


def do_forever(min_loop_time: timedelta) -> Callable:
    """Run the wrapped function in a loop forever.

    The function won't be run more often than once per min_loop_time,
    however if it takes longer to run than min_loop_time, then it will
    be run less often than once per min_loop_time.
    """
    def decorator(function: Callable) -> Callable:
        @wraps(function)
        def wrapper(*args, **kwargs):
            while(True):
                start_time = timezone.now()

                function(*args, **kwargs)

                loop_time = timezone.now() - start_time
                if loop_time < min_loop_time:
                    remaining_time = MIN_LOOP_TIME - loop_time
                    time.sleep(remaining_time.seconds)

        return wrapper
    return decorator


@do_forever(MIN_LOOP_TIME)
def retry_failed_downloader_jobs() -> None:
    """Handle downloader jobs that were marked as a failure."""
    failed_jobs = DownloaderJob.objects.filter(success=False, retried=False)
    handle_downloader_jobs(failed_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_hung_downloader_jobs() -> None:
    """Retry downloader jobs that were started but never finished."""
    minimum_start_time = timezone.now() - MAX_DOWNLOADER_RUN_TIME
    potentially_hung_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__lt=minimum_start_time
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=5)
    hung_jobs = []
    for job in potentially_hung_jobs:
        try:
            job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            if job_status != "running":
                # Make sure it didn't finish since our original query.
                job.refresh_from_db()
                if job.end_time is None:
                    hung_jobs.append(job)
        except URLNotFoundNomadException:
            hung_jobs.append(job)

    handle_downloader_jobs(hung_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_lost_downloader_jobs() -> None:
    """Retry downloader jobs that went too long without being started.

    Idea: at some point this function could integrate with the spot
    instances to determine if jobs are hanging due to a lack of
    instances. A naive time-based implementation like this could end
    up retrying every single queued job if there were a long period
    during which the price of spot instance is higher than our bid
    price.
    """
    minimum_creation_time = timezone.now() - MAX_QUEUE_TIME
    potentially_lost_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        created_at__lt=minimum_creation_time
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=5)
    lost_jobs = []
    for job in potentially_lost_jobs:
        try:
            job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            # If the job is still pending, then it makes sense that it hasn't started.
            if job_status is not "pending":
                # However if it's not pending, then it may have
                # started since our original query.
                job.refresh_from_db()
                if job.start_time is None:
                    # Nope, this job is lost.
                    lost_jobs.append(job)
        except URLNotFoundNomadException:
            lost_jobs.append(job)

    handle_downloader_jobs(lost_jobs)


@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    new_job = ProcessorJob(num_retries=num_retries,
                            pipeline_applied=last_job.pipeline_applied)
    new_job.save()

    for original_file in last_job.original_files.all():
        ProcessorJobOriginalFileAssociation.objects.get_or_create(processor_job=new_job,
                                                          original_file=original_file)

    for data_set in last_job.data_sets.all():
        ProcessorJobDataSetAssociation.objects.get_or_create(processor_job=new_job,
                                                     data_set=data_set)

    logger.info("Requeuing Processor Job which had ID %d with a new Processor Job with ID %d.",
                last_job.id,
                new_job.id)
    send_job(ProcessorPipeline[last_job.pipeline_applied], new_job)

    last_job.retried = True
    last_job.success = False
    last_job.retried_job = new_job
    last_job.save()


def handle_processor_jobs(jobs: List[ProcessorJob]) -> None:
    """For each job in jobs, either retry it or log it."""
    for job in jobs:
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_processor_job(job)
        else:
            handle_repeated_failure(job)


@do_forever(MIN_LOOP_TIME)
def retry_failed_processor_jobs() -> None:
    """Handle processor jobs that were marked as a failure."""
    failed_jobs = ProcessorJob.objects.filter(success=False, retried=False)
    handle_processor_jobs(failed_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_hung_processor_jobs() -> None:
    """Retry processor jobs that were started but never finished."""
    minimum_start_time = timezone.now() - MAX_PROCESSOR_RUN_TIME
    potentially_hung_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__lt=minimum_start_time
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=5)
    hung_jobs = []
    for job in potentially_hung_jobs:
        try:
            job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            if job_status != "running":
                # Make sure it didn't finish since our original query.
                job.refresh_from_db()
                if job.end_time is None:
                    hung_jobs.append(job)
        except URLNotFoundNomadException:
            hung_jobs.append(job)

    handle_processor_jobs(hung_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_lost_processor_jobs() -> None:
    """Retry processor jobs which never even got started for too long."""
    minimum_creation_time = timezone.now() - MAX_QUEUE_TIME
    potentially_lost_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        created_at__lt=minimum_creation_time
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=5)
    lost_jobs = []
    for job in potentially_lost_jobs:
        try:
            job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            # If the job is still pending, then it makes sense that it hasn't started.
            if job_status is not "pending":
                # However if it's not pending, then it may have
                # started since our original query.
                job.refresh_from_db()
                if job.start_time is None:
                    # Nope, this job is lost.
                    lost_jobs.append(job)
        except URLNotFoundNomadException:
            lost_jobs.append(job)

    handle_processor_jobs(lost_jobs)


def monitor_jobs():
    """Runs a thread for each job monitoring loop."""
    functions = [retry_failed_downloader_jobs,
                 retry_hung_downloader_jobs,
                 retry_lost_downloader_jobs,
                 retry_failed_processor_jobs,
                 retry_hung_processor_jobs,
                 retry_lost_processor_jobs]

    threads = []
    for f in functions:
        thread = Thread(target=f, name=f.__name__)
        thread.start()
        threads.append(thread)

    # Make sure that no threads die quietly.
    while(True):
        for thread in threads:
            thread.join(THREAD_WAIT_TIME)
            if not thread.is_alive():
                logger.error("Foreman Thread for the function %s has died!!!!", thread.name)

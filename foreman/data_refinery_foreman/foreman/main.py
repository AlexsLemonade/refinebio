import time
from typing import Callable, List
from threading import Thread
from functools import wraps
from retrying import retry
from datetime import timedelta
from django.utils import timezone
from django.db import transaction
from data_refinery_common.models import (
    WorkerJob,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

# For now it seems like no jobs should take longer than a day to be
# either picked up or run.
MAX_DOWNLOADER_RUN_TIME = timedelta(days=1)
MAX_PROCESSOR_RUN_TIME = timedelta(days=1)
MAX_QUEUE_TIME = timedelta(days=1)

# To prevent excessive spinning loop no more than once every 10
# seconds.
MIN_LOOP_TIME = timedelta(seconds=10)
THREAD_WAIT_TIME = 10.0


@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_downloader_job(last_job: DownloaderJob) -> None:
    """Queues a new downloader job.

    The new downloader job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    new_job = DownloaderJob.create_job_and_relationships(num_retries=num_retries,
                                                         batches=list(last_job.batches.all()),
                                                         downloader_task=last_job.downloader_task)
    logger.info("Requeuing Downloader Job which had ID %d with a new Downloader Job with ID %d.",
                last_job.id,
                new_job.id)
    send_job(Downloaders[last_job.downloader_task], new_job.id)

    last_job.retried = True
    last_job.success = False
    last_job.retried_job = new_job
    last_job.save()


def handle_repeated_failure(job: WorkerJob) -> None:
    """If a job fails too many times, log it and stop retrying."""
    # Not strictly retried but will prevent the job from getting
    # retried any more times.
    job.retried = True

    # success may already be False, but if it was a hung or lost job
    # this will ensure it's marked as failed.
    job.success = False
    job.save()

    # At some point this should become more noisy/attention
    # grabbing. However for the time just logging should be sufficient
    # because all log messages will be closely monitored during early
    # testing stages.
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
    hung_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__lt=minimum_start_time
    )

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
    lost_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        created_at__lt=minimum_creation_time
    )

    handle_downloader_jobs(lost_jobs)


@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    new_job = ProcessorJob.create_job_and_relationships(num_retries=num_retries,
                                                        batches=list(last_job.batches.all()),
                                                        pipeline_applied=last_job.pipeline_applied)
    logger.info("Requeuing Processor Job which had ID %d with a new Processor Job with ID %d.",
                last_job.id,
                new_job.id)
    send_job(ProcessorPipeline[last_job.pipeline_applied], new_job.id)

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
    hung_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__lt=minimum_start_time
    )

    handle_processor_jobs(hung_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_lost_processor_jobs() -> None:
    """Retry processor jobs who never even got started for too long."""
    minimum_creation_time = timezone.now() - MAX_QUEUE_TIME
    lost_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        created_at__lt=minimum_creation_time
        # TEMPORARY for Jackie's grant
    ).exclude(pipeline_applied=ProcessorPipeline.NONE.value)

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

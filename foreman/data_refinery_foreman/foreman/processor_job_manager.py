from typing import List

import data_refinery_foreman.foreman.utils as utils
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import get_capacity_for_jobs, send_job
from data_refinery_common.models import ProcessorJob
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_foreman.foreman.job_requeuing import requeue_processor_job

logger = get_and_configure_logger(__name__)


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

        if job.num_retries < utils.MAX_NUM_RETRIES:
            if requeue_processor_job(job):
                jobs_dispatched = jobs_dispatched + 1
        else:
            utils.handle_repeated_failure(job)


def retry_failed_processor_jobs() -> None:
    """Handle processor jobs that were marked as a failure.

    Ignores Janitor jobs since they are queued every half hour anyway."""
    failed_jobs = (
        ProcessorJob.failed_objects.filter(created_at__gt=utils.JOB_CREATED_AT_CUTOFF)
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

    if queue_capacity <= 0:
        logger.info(
            "Not handling failed (explicitly-marked-as-failure) processor jobs "
            "because there is no capacity for them."
        )

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
        ProcessorJob.hung_objects.filter(created_at__gt=utils.JOB_CREATED_AT_CUTOFF)
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )
    paginator = Paginator(potentially_hung_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling hung (started-but-never-finished) processor jobs "
            "because there is no capacity for them."
        )

    while queue_capacity > 0:
        hung_jobs = utils.check_hung_jobs(database_page.object_list)

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


def retry_lost_processor_jobs() -> None:
    """Retry processor jobs that were started but never finished."""
    potentially_lost_jobs = (
        ProcessorJob.lost_objects.filter(created_at__gt=utils.JOB_CREATED_AT_CUTOFF)
        .exclude(pipeline_applied="JANITOR")
        .order_by("created_at")
        .prefetch_related("original_files__samples")
    )
    paginator = Paginator(potentially_lost_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling lost (never-started) processor jobs "
            "because there is no capacity for them."
        )

    while queue_capacity > 0:
        lost_jobs = utils.check_lost_jobs(database_page.object_list)

        if lost_jobs:
            logger.info(
                "Handling page %d of lost (never-started) processor jobs!",
                database_page_count,
                jobs_count=len(lost_jobs),
            )
            handle_processor_jobs(lost_jobs)

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = get_capacity_for_jobs()
        else:
            break


def retry_unqueued_processor_jobs() -> None:
    """Retry processor jobs that never made it into the Batch job queue."""
    potentially_lost_jobs = ProcessorJob.unqueued_objects.filter(
        created_at__gt=utils.JOB_CREATED_AT_CUTOFF
    ).order_by("created_at")
    paginator = Paginator(potentially_lost_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling unqueued processor jobs " "because there is no capacity for them."
        )

    while queue_capacity > 0:
        for processor_job in database_page.object_list:
            if send_job(
                ProcessorPipeline[processor_job.pipeline_applied],
                job=processor_job,
                is_dispatch=True,
            ):
                queue_capacity -= 1
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            break

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = get_capacity_for_jobs()
        else:
            break

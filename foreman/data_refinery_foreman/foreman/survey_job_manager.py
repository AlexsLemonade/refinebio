from typing import List

import data_refinery_foreman.foreman.utils as utils
from data_refinery_common.job_lookup import SurveyJobTypes
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import SurveyJob
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_foreman.foreman.job_requeuing import requeue_survey_job

logger = get_and_configure_logger(__name__)


def handle_survey_jobs(jobs: List[SurveyJob], queue_capacity: int = None) -> None:
    """For each job in jobs, either retry it or log it.

    No more than queue_capacity jobs will be retried.
    """
    if queue_capacity is None:
        queue_capacity = utils.get_capacity_for_jobs()

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if jobs_dispatched >= queue_capacity:
            logger.info(
                "We hit the maximum total jobs ceiling,"
                " so we're not handling any more survey jobs now."
            )
            return

        if job.num_retries < utils.MAX_NUM_RETRIES:
            requeue_survey_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            utils.handle_repeated_failure(job)


def retry_failed_survey_jobs() -> None:
    """Handle survey jobs that were marked as a failure."""
    failed_jobs = SurveyJob.failed_objects.filter(
        created_at__gt=utils.JOB_CREATED_AT_CUTOFF
    ).order_by("pk")

    paginator = Paginator(failed_jobs, 200)
    page = paginator.page()
    page_count = 0

    if len(page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = utils.get_capacity_for_jobs()

    while queue_capacity > 0:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) survey jobs!", page_count
        )
        handle_survey_jobs(page.object_list, queue_capacity)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
            queue_capacity = utils.get_capacity_for_jobs()
        else:
            break


def retry_hung_survey_jobs() -> None:
    """Retry survey jobs that were started but never finished."""
    potentially_hung_jobs = SurveyJob.hung_objects.filter(
        created_at__gt=utils.JOB_CREATED_AT_CUTOFF
    ).order_by("created_at")
    paginator = Paginator(potentially_hung_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = utils.get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling hung (started-but-never-finished) survey jobs "
            "because there is no capacity for them."
        )
    while queue_capacity > 0:
        hung_jobs = utils.check_hung_jobs(database_page.object_list)

        if hung_jobs:
            logger.info(
                "Handling page %d of hung (started-but-never-finished) survey jobs!",
                database_page_count,
                jobs_count=len(hung_jobs),
            )
            handle_survey_jobs(hung_jobs)

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = utils.get_capacity_for_jobs()
        else:
            break


def retry_lost_survey_jobs() -> None:
    """Retry survey jobs that were started but never finished."""
    potentially_lost_jobs = SurveyJob.lost_objects.filter(
        created_at__gt=utils.JOB_CREATED_AT_CUTOFF
    ).order_by("created_at")
    paginator = Paginator(potentially_lost_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = utils.get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info(
            "Not handling lost (never-started) survey jobs "
            "because there is no capacity for them."
        )
    while queue_capacity > 0:
        lost_jobs = utils.check_lost_jobs(database_page.object_list)

        if lost_jobs:
            logger.info(
                "Handling page %d of lost (never-started) survey jobs!",
                database_page_count,
                jobs_count=len(lost_jobs),
            )
            handle_survey_jobs(lost_jobs)

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = utils.get_capacity_for_jobs()
        else:
            break


def retry_unqueued_survey_jobs() -> None:
    """Retry survey jobs that never made it into the Batch job queue."""
    potentially_lost_jobs = SurveyJob.unqueued_objects.filter(
        created_at__gt=utils.JOB_CREATED_AT_CUTOFF
    ).order_by("created_at")
    paginator = Paginator(potentially_lost_jobs, utils.PAGE_SIZE, "created_at")
    database_page = paginator.page()
    database_page_count = 0

    if len(database_page.object_list) <= 0:
        # No failed jobs, nothing to do!
        return

    queue_capacity = utils.get_capacity_for_jobs()

    if queue_capacity <= 0:
        logger.info("Not handling unqueued survey jobs " "because there is no capacity for them.")

    while queue_capacity > 0:
        for survey_job in database_page.object_list:
            if send_job(SurveyJobTypes.SURVEYOR, job=survey_job, is_dispatch=True):
                queue_capacity -= 1
        else:
            # Can't communicate with Batch just now, leave the job for a later loop.
            break

        if database_page.has_next():
            database_page = paginator.page(database_page.next_page_number())
            database_page_count += 1
            queue_capacity = utils.get_capacity_for_jobs()
        else:
            break

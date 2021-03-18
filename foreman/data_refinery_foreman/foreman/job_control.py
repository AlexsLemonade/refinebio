import datetime
import time

from django.conf import settings
from django.utils import timezone

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import ComputedFile, ProcessorJob
from data_refinery_foreman.foreman.downloader_job_manager import (
    retry_failed_downloader_jobs,
    retry_hung_downloader_jobs,
    retry_lost_downloader_jobs,
)
from data_refinery_foreman.foreman.processor_job_manager import (
    retry_failed_processor_jobs,
    retry_hung_processor_jobs,
    retry_lost_processor_jobs,
)
from data_refinery_foreman.foreman.survey_job_manager import (
    retry_failed_survey_jobs,
    retry_hung_survey_jobs,
    retry_lost_survey_jobs,
)

logger = get_and_configure_logger(__name__)

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


def monitor_jobs():
    """Main Foreman thread that helps manage the Nomad job queue.

    Will find jobs that failed, hung, or got lost and requeue them.

    Also will queue up Janitor jobs regularly to free up disk space.

    Also cleans jobs out of the Nomad queue which cannot be queued
    because the volume containing the job's data isn't mounted.

    It does so on a loop forever that won't spin faster than
    MIN_LOOP_TIME, but it may spin slower than that.
    """
    # last_janitorial_time = timezone.now()
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
            retry_hung_processor_jobs,
            retry_lost_processor_jobs,
            retry_failed_downloader_jobs,
            retry_hung_downloader_jobs,
            retry_lost_downloader_jobs,
            retry_failed_survey_jobs,
            retry_hung_survey_jobs,
            retry_lost_survey_jobs,
        ]

        for function in requeuing_functions_in_order:
            try:
                function()
            except Exception:
                logger.exception("Caught exception in %s: ", function.__name__)

        if settings.RUNNING_IN_CLOUD:
            # Disable this for now because this will trigger regardless of
            # whether or not we have an instance to clean up, which means that
            # we could spin up an instance just to run a janitor job.
            # if timezone.now() - last_janitorial_time > JANITOR_DISPATCH_TIME:
            #     send_janitor_jobs()
            #     last_janitorial_time = timezone.now()

            if timezone.now() - last_dbclean_time > DBCLEAN_TIME:
                clean_database()
                last_dbclean_time = timezone.now()

        loop_time = timezone.now() - start_time
        if loop_time < MIN_LOOP_TIME:
            remaining_time = MIN_LOOP_TIME - loop_time
            if remaining_time.seconds > 0:
                time.sleep(remaining_time.seconds)

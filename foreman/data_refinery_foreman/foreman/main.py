import nomad
import socket
import time
from nomad import Nomad
from nomad.api.exceptions import URLNotFoundNomadException
from typing import Callable, List, Set
from threading import Thread
from functools import wraps
from retrying import retry
from datetime import datetime, timedelta
from django.conf import settings
from django.utils import timezone
from django.db import transaction
from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    SurveyJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation,
    ProcessorJobDatasetAssociation,
    SurveyJobKeyValue
)
from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders, SurveyJobTypes
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable, get_env_variable_gracefully


logger = get_and_configure_logger(__name__)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

# This can be overritten by the env var "MAX_TOTAL_JOBS"
DEFAULT_MAX_JOBS = 20000

# The fastest each thread will repeat its checks.
# Could be slower if the thread takes longer than this to check its jobs.
MIN_LOOP_TIME = timedelta(minutes=2)

# The amount of time the main loop will wait in between checking if
# threads are still alive and then heart beating.
THREAD_WAIT_TIME = timedelta(minutes=10)

# How frequently we dispatch Janitor jobs.
JANITOR_DISPATCH_TIME = timedelta(minutes=30)

# How frequently we clean unplaceable jobs out of the Nomad queue.
CLEANUP_QUEUE_TIME = timedelta(minutes=30)


##
# Utilities
##

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

                try:
                    function(*args, **kwargs)
                except:
                    logger.exception("Exception caught by Foreman while running " + function.__name__)

                loop_time = timezone.now() - start_time
                if loop_time < min_loop_time:
                    remaining_time = MIN_LOOP_TIME - loop_time
                    time.sleep(remaining_time.seconds)

        return wrapper
    return decorator


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


def get_active_volumes() -> Set[str]:
    """Returns a Set of indices for volumes that are currently mounted.

    These can be used to determine which jobs would actually be able
    to be placed if they were queued up.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

    volumes = set()
    for node in nomad_client.nodes.get_nodes():
        node_detail = nomad_client.node.get_node(node["ID"])
        if 'Status' in node_detail and node_detail['Status'] == 'ready' \
           and 'Meta' in node_detail and 'volume_index' in node_detail['Meta']:
            volumes.add(node_detail['Meta']['volume_index'])

    return volumes


##
# Downloaders
##


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

    logger.debug("Requeuing Downloader Job which had ID %d with a new Downloader Job with ID %d.",
                 last_job.id,
                 new_job.id)
    try:
        if send_job(Downloaders[last_job.downloader_task], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
    except:
        logger.error("Failed to requeue Downloader Job which had ID %d with a new Downloader Job with ID %d.",
                     last_job.id,
                     new_job.id)
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()


def handle_downloader_jobs(jobs: List[DownloaderJob]) -> None:
    """For each job in jobs, either retry it or log it."""

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = int(get_env_variable_gracefully("MAX_TOTAL_JOBS", DEFAULT_MAX_JOBS))
    len_all_jobs = len(nomad_client.jobs.get_jobs())
    if len_all_jobs >= MAX_TOTAL_JOBS:
        logger.info("Not requeuing job until we're running fewer jobs.")
        return False

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_downloader_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)

        if (count % 100) == 0:
            len_all_jobs = len(nomad_client.jobs.get_jobs())

        if (jobs_dispatched + len_all_jobs) >= MAX_TOTAL_JOBS:
            logger.info("We hit the maximum total jobs ceiling, so we're not handling any more downloader jobs now.")
            return False

    return True

@do_forever(MIN_LOOP_TIME)
def retry_failed_downloader_jobs() -> None:
    """Handle downloader jobs that were marked as a failure."""
    failed_jobs = DownloaderJob.objects.filter(success=False, retried=False)
    handle_downloader_jobs(failed_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_hung_downloader_jobs() -> None:
    """Retry downloader jobs that were started but never finished."""
    potentially_hung_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__isnull=False,
        no_retry=False
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
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
        except Exception:
            logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

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
    potentially_lost_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        no_retry=False
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    lost_jobs = []
    for job in potentially_lost_jobs:
        try:
            if job.nomad_job_id:
                job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
                # If the job is still pending, then it makes sense that it
                # hasn't started and if it's running then it may not have
                # been able to mark the job record as started yet.
                if job_status != "pending" and job_status != "running":
                    logger.debug(("Determined that a downloader job needs to be requeued because its"
                                  " Nomad Job's status is: %s."),
                                 job_status,
                                 job_id=job.id
                    )
                    lost_jobs.append(job)
            else:
                # If there is no nomad_job_id field set, we could be
                # in the small window where the job was created but
                # hasn't yet gotten a chance to be queued.
                # If this job really should be restarted we'll get it in the next loop.
                if timezone.now() - job.created_at > MIN_LOOP_TIME:
                    lost_jobs.append(job)
        except socket.timeout:
            logger.info("Timeout connecting to Nomad - is Nomad down?", job_id=job.id)
        except URLNotFoundNomadException:
            logger.debug(("Determined that a downloader job needs to be requeued because "
                              "querying for its Nomad job failed: "),
                             job_id=job.id
            )
            lost_jobs.append(job)
        except nomad.api.exceptions.BaseNomadException:
            logger.info("Problem connecting to Nomad - is Nomad down?", job_id=job.id)
        except Exception:
            logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

    handle_downloader_jobs(lost_jobs)


##
# Processors
##


@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    # The Salmon pipeline is quite RAM-sensitive.
    # Try it again with an increased RAM amount, if possible.
    new_ram_amount = last_job.ram_amount
    if last_job.pipeline_applied == "SALMON":
        if new_ram_amount == 4096:
            new_ram_amount = 8192
        elif new_ram_amount == 8192:
            new_ram_amount = 12288
        elif new_ram_amount == 12288:
            new_ram_amount = 16384

    new_job = ProcessorJob(num_retries=num_retries,
                           pipeline_applied=last_job.pipeline_applied,
                           ram_amount=new_ram_amount,
                           volume_index=last_job.volume_index)
    new_job.save()

    for original_file in last_job.original_files.all():
        ProcessorJobOriginalFileAssociation.objects.get_or_create(processor_job=new_job,
                                                          original_file=original_file)

    for dataset in last_job.datasets.all():
        ProcessorJobDatasetAssociation.objects.get_or_create(processor_job=new_job,
                                                     dataset=dataset)

    try:
        logger.debug("Requeuing Processor Job which had ID %d with a new Processor Job with ID %d.",
                     last_job.id,
                     new_job.id)
        if send_job(ProcessorPipeline[last_job.pipeline_applied], job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
    except:
        logger.error("Failed to requeue Processor Job which had ID %d with a new Processor Job with ID %d.",
                     last_job.id,
                     new_job.id)
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()


def handle_processor_jobs(jobs: List[ProcessorJob]) -> None:
    """For each job in jobs, either retry it or log it."""

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = int(get_env_variable_gracefully("MAX_TOTAL_JOBS", DEFAULT_MAX_JOBS))
    len_all_jobs = len(nomad_client.jobs.get_jobs())
    if len_all_jobs >= MAX_TOTAL_JOBS:
        logger.info("Not requeuing job until we're running fewer jobs.")
        return False

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_processor_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)

        if (count % 100) == 0:
            len_all_jobs = len(nomad_client.jobs.get_jobs())

        if (jobs_dispatched + len_all_jobs) >= MAX_TOTAL_JOBS:
            logger.info("We hit the maximum total jobs ceiling, so we're not handling any more processor jobs now.")
            return False

    return True

@do_forever(MIN_LOOP_TIME)
def retry_failed_processor_jobs() -> None:
    """Handle processor jobs that were marked as a failure.

    Ignores Janitor jobs since they are queued every half hour anyway."""
    try:
        active_volumes = get_active_volumes()
    except:
        # If we cannot reach Nomad now then we can wait until a later loop.
        pass

    failed_jobs = ProcessorJob.objects.filter(
        success=False,
        retried=False,
        volume_index__in=active_volumes
    ).exclude(
        pipeline_applied="JANITOR"
    )

    if failed_jobs:
        logger.info(
            "Handling failed (explicitly-marked-as-failure) jobs!",
            jobs=failed_jobs
        )
        handle_processor_jobs(failed_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_hung_processor_jobs() -> None:
    """Retry processor jobs that were started but never finished.

    Ignores Janitor jobs since they are queued every half hour anyway."""
    try:
        active_volumes = get_active_volumes()
    except:
        # If we cannot reach Nomad now then we can wait until a later loop.
        pass

    potentially_hung_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__isnull=False,
        no_retry=False,
        volume_index__in=active_volumes
    ).exclude(
        pipeline_applied="JANITOR"
    )


    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
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
        except TypeError:
            # Almost certainly a python-nomad issue:
            # File "/usr/local/lib/python3.5/dist-packages/nomad/api/job.py", line 63, in get_job
            #   return self.request(id, method="get").json()
            # File "/usr/local/lib/python3.5/dist-packages/nomad/api/base.py", line 74, in request
            #   endpoint = self._endpoint_builder(self.ENDPOINT, *args)
            # File "/usr/local/lib/python3.5/dist-packages/nomad/api/base.py", line 28, in _endpoint_builder
            #   u = "/".join(args)
            # TypeError: sequence item 1: expected str instance, NoneType found
            logger.info("Couldn't query Nomad about Processor Job.", processor_job=job.id)
        except Exception:
            logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

    if hung_jobs:
        logger.info(
            "Handling hung (started-but-never-finished) jobs!",
            jobs=hung_jobs
        )
        handle_processor_jobs(hung_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_lost_processor_jobs() -> None:
    """Retry processor jobs which never even got started for too long.

    Ignores Janitor jobs since they are queued every half hour anyway."""
    try:
        active_volumes = get_active_volumes()
    except:
        # If we cannot reach Nomad now then we can wait until a later loop.
        pass

    potentially_lost_jobs = ProcessorJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        no_retry=False,
        volume_index__in=active_volumes
    ).exclude(
        pipeline_applied="JANITOR"
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=5)
    lost_jobs = []
    for job in potentially_lost_jobs:
        try:
            if job.nomad_job_id:
                job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
                # If the job is still pending, then it makes sense that it
                # hasn't started and if it's running then it may not have
                # been able to mark the job record as started yet.
                if job_status != "pending" and job_status != "running":
                    logger.debug(("Determined that a processor job needs to be requeued because its"
                                  " Nomad Job's status is: %s."),
                                 job_status,
                                 job_id=job.id
                    )
                    lost_jobs.append(job)
            else:
                # If there is no nomad_job_id field set, we could be
                # in the small window where the job was created but
                # hasn't yet gotten a chance to be queued.
                # If this job really should be restarted we'll get it in the next loop.
                if timezone.now() - job.created_at > MIN_LOOP_TIME:
                    lost_jobs.append(job)
        except URLNotFoundNomadException:
            logger.exception(("Determined that a processor job needs to be requeued because "
                              "querying for its Nomad job failed: "),
                             job_id=job.id
            )
            lost_jobs.append(job)
        except Exception:
            logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

    if lost_jobs:
        logger.info(
            "Handling lost (never-started) jobs!",
            jobs=lost_jobs
        )
        handle_processor_jobs(lost_jobs)

##
# Surveyors
##

@retry(stop_max_attempt_number=3)
@transaction.atomic
def requeue_survey_job(last_job: SurveyJob) -> None:
    """Queues a new survey job.

    The new survey job will have num_retries one greater than
    last_job.num_retries.
    """

    lost_jobs = []
    num_retries = last_job.num_retries + 1

    new_job = SurveyJob(num_retries=num_retries,
                        source_type=last_job.source_type
                    )

    if new_job.num_retries == 1:
        new_job.ram_amount = 4096
    elif new_job.num_retries in [2, 3]:
        new_job.ram_amount = 16384
    else:
        new_job.ram_amount = 256

    new_job.save()

    keyvalues = SurveyJobKeyValue.objects.filter(survey_job=last_job)

    for keyvalue in keyvalues:
        SurveyJobKeyValue.objects.get_or_create(survey_job=new_job,
                                                key=keyvalue.key,
                                                value=keyvalue.value,
                                            )

    logger.debug("Requeuing SurveyJob which had ID %d with a new SurveyJob with ID %d.",
                 last_job.id,
                 new_job.id)

    try:
        if send_job(SurveyJobTypes.SURVEYOR, job=new_job, is_dispatch=True):
            last_job.retried = True
            last_job.success = False
            last_job.retried_job = new_job
            last_job.save()
        else:
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
    except:
        logger.error("Failed to requeue Survey Job which had ID %d with a new Surevey Job with ID %d.",
                     last_job.id,
                     new_job.id)
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()

    return True


def handle_survey_jobs(jobs: List[SurveyJob]) -> None:
    """For each job in jobs, either retry it or log it."""

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    # Maximum number of total jobs running at a time.
    # We do this now rather than import time for testing purposes.
    MAX_TOTAL_JOBS = int(get_env_variable_gracefully("MAX_TOTAL_JOBS", DEFAULT_MAX_JOBS))
    len_all_jobs = len(nomad_client.jobs.get_jobs())
    if len_all_jobs >= MAX_TOTAL_JOBS:
        logger.info("Not requeuing job until we're running fewer jobs.")
        return False

    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_survey_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)

        if (count % 100) == 0:
            len_all_jobs = len(nomad_client.jobs.get_jobs())

        if (jobs_dispatched + len_all_jobs) >= MAX_TOTAL_JOBS:
            logger.info("We hit the maximum total jobs ceiling, so we're not handling any more survey jobs now.")
            return False

    return True

@do_forever(MIN_LOOP_TIME)
def retry_failed_survey_jobs() -> None:
    """Handle survey jobs that were marked as a failure."""
    failed_jobs = SurveyJob.objects.filter(success=False, retried=False).order_by('pk')
    if failed_jobs:
        logger.info(
            "Handling failed (explicitly-marked-as-failure) jobs!",
            jobs=failed_jobs
        )
        handle_survey_jobs(failed_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_hung_survey_jobs() -> None:
    """Retry survey jobs that were started but never finished."""
    potentially_hung_jobs = SurveyJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__isnull=False,
        no_retry=False
    ).order_by('pk')

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    hung_jobs = []
    for job in potentially_hung_jobs:
        try:
            # Surveyor jobs didn't always have nomad_job_ids. If they
            # don't have one then by this point they've definitely died.
            if job.nomad_job_id:
                job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            else:
                job_status = "absent"

            if job_status != "running":
                # Make sure it didn't finish since our original query.
                job.refresh_from_db()
                if job.end_time is None:
                    hung_jobs.append(job)
        except URLNotFoundNomadException:
            hung_jobs.append(job)
        except Exception:
            logger.exception("Couldn't query Nomad about SurveyJob Job.", survey_job=job.id)

    if hung_jobs:
        logger.info(
            "Handling hung (started-but-never-finished) jobs!",
            jobs=hung_jobs
        )
        handle_survey_jobs(hung_jobs)


@do_forever(MIN_LOOP_TIME)
def retry_lost_survey_jobs() -> None:
    """Retry survey jobs which never even got started for too long."""
    potentially_lost_jobs = SurveyJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        no_retry=False
    ).order_by('pk')

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
    lost_jobs = []

    for job in potentially_lost_jobs:
        try:
            # Surveyor jobs didn't always have nomad_job_ids. If they
            # don't have one then by this point they've definitely died.
            if job.nomad_job_id:
                job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
            else:
                job_status = "absent"

            # If the job is still pending, then it makes sense that it
            # hasn't started and if it's running then it may not have
            # been able to mark the job record as started yet.
            if job_status != "pending" and job_status != "running":
                logger.debug(("Determined that a survey job needs to be requeued because its"
                             " Nomad Job's status is: %s."),
                            job_status,
                            job_id=job.id
                )
                lost_jobs.append(job)
        except URLNotFoundNomadException:
            logger.debug(("Determined that a survey job needs to be requeued because "
                          "querying for its Nomad job failed."),
                         job_id=job.id
            )
            lost_jobs.append(job)
        except Exception:
            logger.exception("Couldn't query Nomad about Processor Job.", survey_job=job.id)

    if lost_jobs:
        logger.info(
            "Handling lost (never-started) jobs!",
            jobs=lost_jobs
        )
        handle_survey_jobs(lost_jobs)

##
# Janitor
##

@do_forever(JANITOR_DISPATCH_TIME)
def send_janitor_jobs():
    """Dispatch a Janitor job for each instance in the cluster"""
    try:
        active_volumes = get_active_volumes()
    except:
        # If we cannot reach Nomad now then we can wait until a later loop.
        pass

    for volume_index in active_volumes:
        new_job = ProcessorJob(num_retries=0,
                               pipeline_applied="JANITOR",
                               ram_amount=2048,
                               volume_index=volume_index)
        new_job.save()
        logger.info("Sending Janitor with index: ",
            job_id=new_job.id,
            index=volume_index
        )
        try:
            send_job(ProcessorPipeline["JANITOR"], job=new_job, is_dispatch=True)
        except Exception as e:
            # If we can't dispatch this job, something else has gone wrong.
            continue


##
# Handling of node cycling
##

@do_forever(CLEANUP_QUEUE_TIME)
def cleanup_the_queue():
    """This cleans up any jobs which cannot currently be queued.

    We often have more volumes than instances because we have enough
    volumes for the scenario where the entire cluster is using the
    smallest instance type, however that doesn't happen very
    often. Therefore it's possible for some volumes to not be mounted,
    which means that jobs which are constrained to run on instances
    with those volumes cannot be placed and just clog up the queue.

    Therefore we clear out jobs of that type every once in a while so
    our queue is dedicated to jobs that can actually be placed.
    """
    # Smasher and QN Reference jobs aren't tied to a specific EBS volume.
    indexed_job_types = [e.value for e in ProcessorPipeline if e.value not in ["SMASHER", "QN_REFERENCE"]]

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

    try:
        active_volumes = get_active_volumes()
        jobs = nomad_client.jobs.get_jobs()
    except:
        # If we cannot reach Nomad now then we can wait until a later loop.
        pass


    jobs_to_kill = []
    for job in jobs:
        # Skip over the Parameterized Jobs because we need those to
        # always be running.
        if "ParameterizedJob" not in job or not job["ParameterizedJob"]:
            continue

        for job_type in indexed_job_types:
            # We're only concerned with jobs that have to be tied to a volume index.
            if "ParentID" not in job or not job["ParentID"].startswith(job_type):
                continue

            # If this job has an index, then its ParentID will
            # have the pattern of <job-type>_<index>_<RAM-amount>
            # and we want to check the value of <index>:
            split_parent_id = job["ParentID"].split("_")
            if len(split_parent_id) < 2:
                continue
            else:
                index = split_parent_id[-2]

            if index not in active_volumes:
                # The index for this job isn't currently mounted, kill
                # the job and decrement the retry counter (since it
                # will be incremented when it is requeued).
                try:
                    nomad_client.job.deregister_job(job["ID"], purge=True)
                    job_record = ProcessorJob(nomad_job_id=job["ID"])
                    job_record.num_retries = job_record.num_retries - 1
                    job_record.save()
                except:
                    # If we can't do this for some reason, we'll get it next loop.
                    pass



##
# Main loop
##

def monitor_jobs():
    """Runs a thread for each job monitoring loop."""

    threads = []

    # Start the thread to dispatch Janitor jobs.
    thread = Thread(target=send_janitor_jobs, name="send_janitor_jobs")
    thread.start()
    threads.append(thread)
    logger.info("Thread started for function: send_janitor_jobs")

    # Start the thread to keep the Nomad clean of jobs it cannot place.
    thread = Thread(target=cleanup_the_queue, name="cleanup_the_queue")
    thread.start()
    threads.append(thread)
    logger.info("Thread started for function: cleanup_the_queue")

    processor_functions = [ retry_failed_processor_jobs,
                            retry_hung_processor_jobs,
                            retry_lost_processor_jobs]

    for f in processor_functions:
        thread = Thread(target=f, name=f.__name__)
        thread.start()
        threads.append(thread)
        logger.info("Thread started for monitoring function: %s", f.__name__)

    # This is only a concern when running at scale.
    if settings.RUNNING_IN_CLOUD:
        # We start the processor threads first so that we don't
        # accidentally queue too many downloader jobs and knock down our
        # source databases. They may take a while to run, and this
        # function only runs once per deploy, so give a generous amount of
        # time, say 5 minutes:
        time.sleep(60*5)

    downloader_functions = [retry_failed_downloader_jobs,
                            retry_hung_downloader_jobs,
                            retry_lost_downloader_jobs]

    for f in downloader_functions:
        thread = Thread(target=f, name=f.__name__)
        thread.start()
        threads.append(thread)
        logger.info("Thread started for monitoring function: %s", f.__name__)

    survey_functions = [retry_failed_survey_jobs,
                            retry_hung_survey_jobs,
                            retry_lost_survey_jobs]

    for f in survey_functions:
        thread = Thread(target=f, name=f.__name__)
        thread.start()
        threads.append(thread)
        logger.info("Thread started for monitoring function: %s", f.__name__)

    # Make sure that no threads die quietly.
    while(True):

        # Write the health file for Monit to check
        now_secs = int(time.time())
        with open('/tmp/foreman_last_time', 'w') as timefile:
            timefile.write(str(now_secs))

        start_time = timezone.now()

        for thread in threads:
            if not thread.is_alive():
                logger.error("Foreman Thread for the function %s has died!!!!", thread.name)

        loop_time = timezone.now() - start_time
        if loop_time < THREAD_WAIT_TIME:
            remaining_time = THREAD_WAIT_TIME - loop_time
            time.sleep(remaining_time.seconds)

        logger.info("The Foreman's heart is beating, but he does not feel.")

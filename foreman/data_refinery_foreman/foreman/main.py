import datetime
import nomad
import socket
import time
import traceback

from django.conf import settings
#from django.core.paginator import Paginator
from data_refinery_foreman.foreman.performant_pagination.pagination import PerformantPaginator as Paginator
from django.db import transaction
from django.utils import timezone
from functools import wraps
from nomad import Nomad
from nomad.api.exceptions import URLNotFoundNomadException
from typing import List, Set

from data_refinery_common.job_lookup import (
    Downloaders,
    ProcessorPipeline,
    SurveyJobTypes,
    does_processor_job_have_samples,
    is_file_rnaseq,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    SurveyJobKeyValue
)
from data_refinery_common.utils import (
    get_active_volumes,
    get_env_variable,
    get_env_variable_gracefully
)


logger = get_and_configure_logger(__name__)

# Maximum number of retries, so the number of attempts will be one
# greater than this because of the first attempt
MAX_NUM_RETRIES = 2

# This can be overritten by the env var "MAX_TOTAL_JOBS"
DEFAULT_MAX_JOBS = 20000

# This is the maximum number of non-dead nomad jobs that can be in the queue.
MAX_TOTAL_DOWNLOADER_JOBS = 0
DOWNLOADER_JOBS_PER_NODE = 40
TIME_OF_LAST_SIZE_CHECK = timezone.now()

# The minimum amount of time in between each iteration of the main
# loop. We could loop much less frequently than every two minutes if
# the work we do takes longer than 2 minutes, but this will prevent
# excessive spinning.
MIN_LOOP_TIME = datetime.timedelta(minutes=2)

# How frequently we dispatch Janitor jobs and clean unplaceable jobs
# out of the Nomad queue.
JANITOR_DISPATCH_TIME = datetime.timedelta(minutes=30)

# This time is currently set so far in the past that it's not doing
# anything. However, should we ever want to suspend work on the
# whatever we already have in our queues/database (perhaps to be able
# to force "important work" to get done), setting this to a recent
# date will prevent the Foreman from queuing/requeuing jobs created
# before this cutoff.
JOB_CREATED_AT_CUTOFF = datetime.datetime(2017, 2, 5, tzinfo=timezone.utc)


def read_config_list(config_file: str) -> List[str]:
    """
    Reads a file and returns a list with one item per line.
    """
    return_list = []
    with open(config_file) as config_list:
        for line in config_list:
            return_list.append(line.strip())

    return return_list


# These are lists of accessions that are for
# subsets of experiments that we want to prioritize.
PEDIATRIC_ACCESSION_LIST = read_config_list("config/pediatric_accessions.txt")
HGU133PLUS2_ACCESSION_LIST = read_config_list("config/hgu133plus2_accessions.txt")


##
# Utilities
##


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

def get_max_downloader_jobs(window=datetime.timedelta(minutes=2), nomad_client=None):
    """
    Fetches the desired maximum number of downloader jobs available based on the cluster size.
    If this has been calculated recently, returns a cached value, else it will calculate it fresh every `window`.
    """
    global MAX_TOTAL_DOWNLOADER_JOBS
    global TIME_OF_LAST_SIZE_CHECK

    if (timezone.now() - TIME_OF_LAST_SIZE_CHECK > window):
        # Assuming they're all similar to an X1, give 50 downloaders per node.
        if not nomad_client:
            nomad_host = get_env_variable("NOMAD_HOST")
            nomad_port = get_env_variable("NOMAD_PORT", "4646")
            nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

        # Minus one because the smasher doesn't run DLs
        MAX_TOTAL_DOWNLOADER_JOBS = (len(nomad_client.nodes) - 1) * DOWNLOADER_JOBS_PER_NODE
        TIME_OF_LAST_SIZE_CHECK = timezone.now()

    if MAX_TOTAL_DOWNLOADER_JOBS > 1000:
        return 1000
    else:
        return MAX_TOTAL_DOWNLOADER_JOBS

##
# Job Prioritization
##


def prioritize_salmon_jobs(jobs: List) -> List:
    """Prioritizes salmon experiments based on how close to completion they are.

    This is because salmon experiments have a final processing step
    that must be performed on all the samples in the experiment, so if
    9/10 samples in an experiment are processed then they can't
    actually be used until that last sample is processed.
    """
    # The strategy for doing so is to build a mapping between every
    # salmon job in `jobs` to a priority. This priority will be what
    # percentage of the samples in this experiment have been
    # processed. Once we have that mapping we can sort those jobs by
    # that priority and move them to the front of the list.
    prioritized_jobs = []
    for job in jobs:
        try:
            if type(job) == ProcessorJob and not does_processor_job_have_samples(job):
                continue

            # Salmon jobs are specifc to one sample.
            sample = job.get_samples().pop()

            # Skip jobs that aren't for Salmon. Handle both ProcessorJobs and DownloaderJobs.
            if type(job) is ProcessorJob and job.pipeline_applied != ProcessorPipeline.SALMON.value:
                continue
            elif type(job) is DownloaderJob:
                is_salmon_sample = False
                for original_file in sample.original_files.all():
                    if is_file_rnaseq(original_file.filename):
                        is_salmon_sample = True

                if not is_salmon_sample:
                    continue

            # Get a set of unique samples that share at least one
            # experiment with the sample this job is for.
            related_samples = set()
            for experiment in sample.experiments.all():
                for related_sample in experiment.samples.all():
                    related_samples.add(related_sample)

            # We cannot simply filter on is_processed because that field
            # doesn't get set until every sample in an experiment is processed.
            # Instead we are looking for one successful processor job.
            processed_samples = 0
            for related_sample in related_samples:
                original_files = related_sample.original_files
                if original_files.count() == 0:
                    logger.error("Salmon sample found without any original files!!!", sample=related_sample)
                elif original_files.first().processor_jobs.filter(success=True).count() >= 1:
                    processed_samples += 1

            experiment_completion_percent = processed_samples / len(related_samples)
            prioritized_jobs.append({"job": job, "priority": experiment_completion_percent})
        except:
            logger.debug("Exception caught while prioritizing salmon jobs!", job=job)


    sorted_job_mappings = sorted(prioritized_jobs, reverse=True, key=lambda k: k["priority"])
    sorted_jobs = [job_mapping["job"] for job_mapping in sorted_job_mappings]

    # Remove all the jobs we're moving to the front of the list
    for job in sorted_jobs:
        jobs.remove(job)

    return sorted_jobs + jobs


def prioritize_zebrafish_jobs(jobs: List) -> List:
    """Moves zebrafish jobs to the beginnging of the input list."""
    zebrafish_jobs = []
    for job in jobs:
        try:
            if type(job) == ProcessorJob and not does_processor_job_have_samples(job):
                continue

            # There aren't cross-species jobs, so just checking one sample's organism will be sufficient.
            samples = job.get_samples()

            for sample in samples:
                if sample.organism.name == 'DANIO_RERIO':
                    zebrafish_jobs.append(job)
                    break
        except:
            logger.debug("Exception caught while prioritizing zebrafish jobs!", job=job)

    # Remove all the jobs we're moving to the front of the list
    for job in zebrafish_jobs:
        jobs.remove(job)

    return zebrafish_jobs + jobs


def prioritize_jobs_by_accession(jobs: List, accession_list: List[str]) -> List:
    """Moves jobs whose accessions are in accession_lst to the beginning of the input list."""
    prioritized_jobs = []
    for job in jobs:
        try:
            if type(job) == ProcessorJob and not does_processor_job_have_samples(job):
                continue

            # All samples in a job correspond to the same experiment, so just check one sample.
            samples = job.get_samples()

            # Iterate through all the samples' experiments until one is
            # found with an accession in `accession_list`.
            is_prioritized_job = False
            for sample in samples:
                if is_prioritized_job:
                    # We found one! So stop looping
                    break

                for experiment in sample.experiments.all():
                    if experiment.accession_code in accession_list:
                        prioritized_jobs.append(job)
                        is_prioritized_job = True
                        break
        except:
            logger.exception("Exception caught while prioritizing jobs by accession!", job=job)

    # Remove all the jobs we're moving to the front of the list
    for job in prioritized_jobs:
        jobs.remove(job)

    return prioritized_jobs + jobs


##
# Downloaders
##

def requeue_downloader_job(last_job: DownloaderJob) -> None:
    """Queues a new downloader job.

    The new downloader job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    ram_amount = last_job.ram_amount
    if last_job.failure_reason is not None and 'harakiri' not in last_job.failure_reason:
        if ram_amount == 1024:
            ram_amount = 4096
        elif ram_amount == 4096:
            ram_amount = 8192

    new_job = DownloaderJob(num_retries=num_retries,
                            downloader_task=last_job.downloader_task,
                            ram_amount=ram_amount,
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


def count_downloader_jobs_in_queue(nomad_client: Nomad) -> int:
    """Counts how many downloader jobs in the Nomad queue do not have status of 'dead'."""
    all_downloader_jobs = nomad_client.jobs.get_jobs(prefix="DOWNLOADER")

    total = 0
    for job in all_downloader_jobs:
        if job['ParameterizedJob']:
            total = total + job['JobSummary']['Children']['Pending']
            total = total + job['JobSummary']['Children']['Running']

    return total


def handle_downloader_jobs(jobs: List[DownloaderJob]) -> None:
    """For each job in jobs, either retry it or log it."""

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    num_downloader_jobs = count_downloader_jobs_in_queue(nomad_client)
    current_max_downloader_jobs = get_max_downloader_jobs(nomad_client=nomad_client)
    if num_downloader_jobs >= current_max_downloader_jobs:
        logger.info("Not requeuing job until we're running fewer jobs.")
        return False

    # We want zebrafish data first, then hgu133plus2, then data
    # related to pediatric cancer, then to finish salmon experiments
    # that are close to completion.
    # Each function moves the jobs it prioritizes to the front of the
    # list, so apply them in backwards order.
    # jobs = prioritize_salmon_jobs(jobs)
    # jobs = prioritize_jobs_by_accession(jobs, PEDIATRIC_ACCESSION_LIST)
    # jobs = prioritize_jobs_by_accession(jobs, HGU133PLUS2_ACCESSION_LIST)
    # jobs = prioritize_zebrafish_jobs(jobs)

    num_to_dispatch = current_max_downloader_jobs - num_downloader_jobs
    jobs_dispatched = 0
    for count, job in enumerate(jobs):
        if job.num_retries < MAX_NUM_RETRIES:
            requeue_downloader_job(job)
            jobs_dispatched = jobs_dispatched + 1
        else:
            handle_repeated_failure(job)

        if jobs_dispatched > num_to_dispatch:
            break

    return True

def retry_failed_downloader_jobs() -> None:
    """Handle downloader jobs that were marked as a failure."""
    failed_jobs = DownloaderJob.objects.filter(
        success=False,
        retried=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )

    if failed_jobs.count() == 0:
        return

    paginator = Paginator(failed_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) downloader jobs!",
            page_count

        )
        handle_downloader_jobs(page.object_list)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

def retry_hung_downloader_jobs() -> None:
    """Retry downloader jobs that were started but never finished."""
    potentially_hung_jobs = DownloaderJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__isnull=False,
        no_retry=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    if potentially_hung_jobs.count() == 0:
        return

    paginator = Paginator(potentially_hung_jobs, 200)
    page = paginator.page()

    page_count = 0
    while True:
        hung_jobs = []
        for job in page.object_list:
            try:
                job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]
                if job_status != "running":
                    # Make sure it didn't finish since our original query.
                    job.refresh_from_db()
                    if job.end_time is None:
                        hung_jobs.append(job)
            except URLNotFoundNomadException:
                hung_jobs.append(job)
            except nomad.api.exceptions.BaseNomadException:
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about Downloader Job.", downloader_job=job.id)

        if hung_jobs:
            logger.info(
                "Handling page %d of hung (started-but-never-finished) downloader jobs!",
                page_count,
                jobs_count=len(hung_jobs)
            )
            handle_downloader_jobs(hung_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

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
        no_retry=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    if potentially_lost_jobs.count() == 0:
        return

    paginator = Paginator(potentially_lost_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        lost_jobs = []
        for job in page.object_list:
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
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about Downloader Job.", downloader_job=job.id)

        if lost_jobs:
            logger.info(
                "Handling page %d of lost (never-started) downloader jobs!",
                page_count,
                len_jobs=len(lost_jobs)
            )
            handle_downloader_jobs(lost_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break


##
# Processors
##

def requeue_processor_job(last_job: ProcessorJob) -> None:
    """Queues a new processor job.

    The new processor job will have num_retries one greater than
    last_job.num_retries.
    """
    num_retries = last_job.num_retries + 1

    # The Salmon pipeline is quite RAM-sensitive.
    # Try it again with an increased RAM amount, if possible.
    new_ram_amount = last_job.ram_amount

    # These initial values are set in common/job_lookup.py:determine_ram_amount
    if last_job.pipeline_applied == "SALMON":
        if new_ram_amount == 12288:
            new_ram_amount = 16384
        elif new_ram_amount == 16384:
            new_ram_amount = 32768
        elif new_ram_amount == 32768:
            new_ram_amount = 65536
    # The AFFY pipeline is somewhat RAM-sensitive.
    # Try it again with an increased RAM amount, if possible.
    elif last_job.pipeline_applied == "AFFY_TO_PCL":
        if new_ram_amount == 2048:
            new_ram_amount = 4096
        elif new_ram_amount == 4096:
            new_ram_amount = 8192

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

    # We want zebrafish data first, then hgu133plus2, then data
    # related to pediatric cancer, then to finish salmon experiments
    # that are close to completion.
    # Each function moves the jobs it prioritizes to the front of the
    # list, so apply them in backwards order.
    # jobs = prioritize_salmon_jobs(jobs)
    # jobs = prioritize_jobs_by_accession(jobs, PEDIATRIC_ACCESSION_LIST)
    # jobs = prioritize_jobs_by_accession(jobs, HGU133PLUS2_ACCESSION_LIST)
    # jobs = prioritize_zebrafish_jobs(jobs)

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
        volume_index__in=active_volumes,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).exclude(
        pipeline_applied="JANITOR"
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )

    if failed_jobs.count() == 0:
        return

    paginator = Paginator(failed_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) processor jobs!",
            page_count
        )
        handle_processor_jobs(page.object_list)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

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
        volume_index__in=active_volumes,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).exclude(
        pipeline_applied="JANITOR"
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )


    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    if potentially_hung_jobs.count() == 0:
        return

    paginator = Paginator(potentially_hung_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        hung_jobs = []
        for job in page.object_list:
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
            except nomad.api.exceptions.BaseNomadException:
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

        if hung_jobs:
            logger.info(
                "Handling hung page %d of (started-but-never-finished) processor jobs!",
                page_count,
                len_jobs=len(hung_jobs)
            )
            handle_processor_jobs(hung_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

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
        volume_index__in=active_volumes,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).exclude(
        pipeline_applied="JANITOR"
    ).order_by(
        'id'
    ).prefetch_related(
        "original_files__samples"
    )

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=5)

    if potentially_lost_jobs.count() == 0:
        return

    paginator = Paginator(potentially_lost_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        lost_jobs = []
        for job in page.object_list:
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
                logger.debug(("Determined that a processor job needs to be requeued because "
                                  "querying for its Nomad job failed: "),
                                 job_id=job.id
                )
                lost_jobs.append(job)
            except nomad.api.exceptions.BaseNomadException:
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about Processor Job.", processor_job=job.id)

        if lost_jobs:
            logger.info(
                "Handling lost page %d of (never-started) processor jobs!",
                page_count,
                len_jobs=len(lost_jobs)
            )
            handle_processor_jobs(lost_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

##
# Surveyors
##

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


def retry_failed_survey_jobs() -> None:
    """Handle survey jobs that were marked as a failure."""
    failed_jobs = SurveyJob.objects.filter(
        success=False,
        retried=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by('pk')

    if failed_jobs.count() == 0:
        return

    paginator = Paginator(failed_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        logger.info(
            "Handling page %d of failed (explicitly-marked-as-failure) survey jobs!",
            page_count
        )
        handle_survey_jobs(page.object_list)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break


def retry_hung_survey_jobs() -> None:
    """Retry survey jobs that were started but never finished."""
    potentially_hung_jobs = SurveyJob.objects.filter(
        success=None,
        retried=False,
        end_time=None,
        start_time__isnull=False,
        no_retry=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by('pk')

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    if potentially_hung_jobs.count() == 0:
        return

    paginator = Paginator(potentially_hung_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        hung_jobs = []
        for job in page.object_list:
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
            except nomad.api.exceptions.BaseNomadException:
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about SurveyJob Job.", survey_job=job.id)

        if hung_jobs:
            logger.info(
                "Handling page %d of hung (started-but-never-finished) survey jobs!",
                page_count,
                len_jobs=len(hung_jobs)
            )
            handle_survey_jobs(hung_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

def retry_lost_survey_jobs() -> None:
    """Retry survey jobs which never even got started for too long."""
    potentially_lost_jobs = SurveyJob.objects.filter(
        success=None,
        retried=False,
        start_time=None,
        end_time=None,
        no_retry=False,
        created_at__gt=JOB_CREATED_AT_CUTOFF
    ).order_by('pk')

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

    if potentially_lost_jobs.count() == 0:
        return

    paginator = Paginator(potentially_lost_jobs, 200)
    page = paginator.page()
    page_count = 0
    while True:
        lost_jobs = []
        for job in page.object_list:
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
            except nomad.api.exceptions.BaseNomadException:
                raise
            except Exception:
                logger.exception("Couldn't query Nomad about Processor Job.", survey_job=job.id)

        if lost_jobs:
            logger.info(
                "Handling page %d of lost (never-started) survey jobs!",
                page_count,
                len_jobs=len(lost_jobs)
            )
            handle_survey_jobs(lost_jobs)

        if page.has_next():
            page = paginator.page(page.next_page_number())
            page_count = page_count + 1
        else:
            break

##
# Janitor
##

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
        return


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
    """Main Foreman thread that helps manage the Nomad job queue.

    Will find jobs that failed, hung, or got lost and requeue them.

    Also will queue up Janitor jobs regularly to free up disk space.

    Also cleans jobs out of the Nomad queue which cannot be queued
    because the volume containing the job's data isn't mounted.

    It does so on a loop forever that won't spin faster than
    MIN_LOOP_TIME, but it may spin slower than that.
    """
    last_janitorial_time = None
    while(True):
        # Perform two heartbeats, one for the logs and one for Monit:
        logger.info("The Foreman's heart is beating, but he does not feel.")

        # Write the health file for Monit to check
        now_secs = int(time.time())
        with open('/tmp/foreman_last_time', 'w') as timefile:
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
            retry_lost_survey_jobs
        ]

        for function in requeuing_functions_in_order:
            try:
                function()
            except Exception as e:
                logger.error("Caught exception in %s: ", function.__name__)
                traceback.print_exc(chain=False)

        if not last_janitorial_time \
           or timezone.now() - last_janitorial_time > JANITOR_DISPATCH_TIME:
            send_janitor_jobs()
            cleanup_the_queue()
            last_janitorial_time = timezone.now()

        loop_time = timezone.now() - start_time
        if loop_time < MIN_LOOP_TIME:
            remaining_time = MIN_LOOP_TIME - loop_time
            if remaining_time.seconds > 0:
                time.sleep(remaining_time.seconds)

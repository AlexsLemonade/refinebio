import datetime
import psutil
import signal
import sys

from django.db import transaction
from django.conf import settings
from django.utils import timezone
from retrying import retry
from typing import List, Dict

from data_refinery_common.job_lookup import ProcessorPipeline, determine_processor_pipeline, determine_ram_amount
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
)
from data_refinery_common.utils import (
    get_env_variable,
    get_instance_id,
)

logger = get_and_configure_logger(__name__)
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")
# TODO: extend this list.
BLACKLISTED_EXTENSIONS = ["xml", "chp", "exp"]
CURRENT_JOB = None

def get_max_jobs_for_current_node():
    """ Determine the maximum number of Downloader jobs that this node should sustain,
    based on total system RAM."""

    total_vm = psutil.virtual_memory().total
    gb = int(total_vm / 1000000000)
    logger.info("Detected " + str(gb) + "GB of RAM.")

    # We basically want to hit 2GB/s total across our entire cluster. Each job hits ~20MB/s,
    # so it'd take 100 jobs to hit our limit. Our cluster has 10TB of RAM, which is ~10000 GB.
    # 100 jobs per 10000 GB = 1 job per 100 GB
    max_jobs = (gb/100)

    # However this will make sure we can still run a few jobs in local environments and in CI.
    if max_jobs < 5:
        max_jobs = 5

    return max_jobs

MAX_DOWNLOADER_JOBS_PER_NODE = get_max_jobs_for_current_node()

def signal_handler(sig, frame):
    """Signal Handler, works for both SIGTERM and SIGINT"""
    global CURRENT_JOB
    if not CURRENT_JOB:
        sys.exit(0)
    else:
        CURRENT_JOB.start_time = None
        CURRENT_JOB.num_retries = CURRENT_JOB.num_retries - 1
        CURRENT_JOB.failure_reason = "Interruped by SIGTERM/SIGINT"
        CURRENT_JOB.save()
        sys.exit(0)

def start_job(job_id: int, max_downloader_jobs_per_node=MAX_DOWNLOADER_JOBS_PER_NODE, force_harakiri=False) -> DownloaderJob:
    """Record in the database that this job is being started.

    Retrieves the job from the database and returns it after marking
    it as started.
    """
    logger.debug("Starting Downloader Job.", downloader_job=job_id)
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except DownloaderJob.DoesNotExist:
        logger.error("Cannot find downloader job record.", downloader_job=job_id)
        raise

    worker_id = get_instance_id()
    num_downloader_jobs_currently_running = DownloaderJob.objects.filter(
                                worker_id=worker_id,
                                start_time__isnull=False,
                                end_time__isnull=True,
                                success__isnull=True,
                                retried=False
                            ).count()

    # Death and rebirth.
    # if settings.RUNNING_IN_CLOUD or force_harakiri:
    #     if num_downloader_jobs_currently_running >= int(max_downloader_jobs_per_node):
    #         # Wait for the death window
    #         while True:
    #             seconds = datetime.datetime.now().second
    #             # Mass harakiri happens every 15 seconds.
    #             if seconds % 15 == 0:
    #                 job.start_time = None
    #                 job.num_retries = job.num_retries - 1
    #                 job.failure_reason = "Killed by harakiri"
    #                 job.success = False
    #                 job.save()

    #                 # What is dead may never die!
    #                 sys.exit(0)

    # This job should not have been started.
    if job.start_time is not None:
        logger.error("This downloader job has already been started!!!", downloader_job=job.id)
        raise Exception("downloaders.start_job called on a job that has already been started!")

    original_file = job.original_files.first()
    if original_file and original_file.has_been_processed():
        logger.error(("Sample has a good computed file, it must have been processed, "
                      "so it doesn't need to be downloaded! Aborting!"),
                     job_id=job.id,
                     original_file_id=original_file.id
        )
        job.start_time = timezone.now()
        job.failure_reason = "Was told to redownload a successfully proccessed file."
        job.success = False
        job.no_retry = True
        job.end_time = timezone.now()
        job.save()
        sys.exit(0)

    # Set up the SIGTERM handler so we can appropriately handle being interrupted.
    # (`docker stop` uses SIGTERM, not SIGINT.)
    # (however, Nomad sends an SIGINT so catch both.)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    job.worker_id = worker_id
    job.worker_version = SYSTEM_VERSION
    job.start_time = timezone.now()
    job.save()

    global CURRENT_JOB
    CURRENT_JOB = job

    return job


def end_downloader_job(job: DownloaderJob, success: bool):
    """
    Record in the database that this job has completed.
    """
    if success:
        logger.debug("Downloader Job completed successfully.",
                    downloader_job=job.id)
    else:
        # Should be set by now, but make sure.
        success = False
        file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
        for file_assoc in file_assocs:
            file_assoc.original_file.delete_local_file()
            file_assoc.original_file.is_downloaded = False
            file_assoc.original_file.save()

        if not job.failure_reason:
            logger.error("Downloader job failed without having failure_reason set. FIX ME!!!!!!!!",
                         downloader_job=job.id,
                         downloader_task=job.downloader_task)
        else:
            logger.info("Downloader job failed!",
                        downloader_job=job.id,
                        downloader_task=job.downloader_task,
                        failure_reason=job.failure_reason)

    job.success = success
    job.end_time = timezone.now()
    job.save()

def delete_if_blacklisted(original_file: OriginalFile) -> OriginalFile:
    extension = original_file.filename.split(".")[-1]
    if extension.lower() in BLACKLISTED_EXTENSIONS:
        logger.debug("Original file had a blacklisted extension of %s, skipping",
                     extension,
                     original_file=original_file.id)

        original_file.delete_local_file()
        original_file.is_downloaded = False
        original_file.save()
        return None

    return original_file


def create_processor_jobs_for_original_files(original_files: List[OriginalFile],
                                             downloader_job: DownloaderJob=None):
    """
    Create a processor jobs and queue a processor task for samples related to an experiment.
    """
    for original_file in original_files:
        sample_object = original_file.samples.first()

        if not delete_if_blacklisted(original_file):
            continue

        # Fix for: https://github.com/AlexsLemonade/refinebio/issues/968
        # Basically, we incorrectly detected technology/manufacturers
        # for many Affymetrix samples and this is a good place to fix
        # some of them.
        if original_file.is_affy_data():
            # Only Affymetrix Microarrays produce .CEL files
            sample_object.technology = 'MICROARRAY'
            sample_object.manufacturer = 'AFFYMETRTIX'
            sample_object.save()

        pipeline_to_apply = determine_processor_pipeline(sample_object, original_file)

        if pipeline_to_apply == ProcessorPipeline.NONE:
            logger.info("No valid processor pipeline found to apply to sample.",
                        sample=sample_object.id,
                        original_file=original_files[0].id)
            original_file.delete_local_file()
            original_file.is_downloaded = False
            original_file.save()
        else:
            processor_job = ProcessorJob()
            processor_job.pipeline_applied = pipeline_to_apply.value
            processor_job.ram_amount = determine_ram_amount(sample_object, processor_job)
            processor_job.save()

            assoc = ProcessorJobOriginalFileAssociation()
            assoc.original_file = original_file
            assoc.processor_job = processor_job
            assoc.save()

            if downloader_job:
                logger.debug("Queuing processor job.",
                             processor_job=processor_job.id,
                             original_file=original_file.id,
                             downloader_job=downloader_job.id)
            else:
                logger.debug("Queuing processor job.",
                             processor_job=processor_job.id,
                             original_file=original_file.id)

            try:
                send_job(pipeline_to_apply, processor_job)
            except:
                # If we cannot queue the job now the Foreman will do
                # it later.
                pass


def create_processor_job_for_original_files(original_files: List[OriginalFile],
                                            downloader_job: DownloaderJob=None):
    """
    Create a processor job and queue a processor task for sample related to an experiment.

    """

    # If there's no original files then we've created all the jobs we need to!
    if len(original_files) == 0:
        return

    # For anything that has raw data there should only be one Sample per OriginalFile
    sample_object = original_files[0].samples.first()
    pipeline_to_apply = determine_processor_pipeline(sample_object, original_files[0])

    if pipeline_to_apply == ProcessorPipeline.NONE:
        logger.info("No valid processor pipeline found to apply to sample.",
                    sample=sample_object.id,
                    original_file=original_files[0].id)
        for original_file in original_files:
            original_file.delete_local_file()
            original_file.is_downloaded = False
            original_file.save()
    else:
        processor_job = ProcessorJob()
        processor_job.pipeline_applied = pipeline_to_apply.value
        processor_job.ram_amount = determine_ram_amount(sample_object, processor_job)
        processor_job.save()
        for original_file in original_files:
            assoc = ProcessorJobOriginalFileAssociation()
            assoc.original_file = original_file
            assoc.processor_job = processor_job
            assoc.save()

        logger.debug("Queuing processor job.",
                     processor_job=processor_job.id)

        try:
            send_job(pipeline_to_apply, processor_job)
        except:
            # If we cannot queue the job now the Foreman will do
            # it later.
            pass

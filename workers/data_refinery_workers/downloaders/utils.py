import datetime
import signal
import sys
from typing import Dict, List

from django.conf import settings
from django.db import transaction
from django.utils import timezone

import psutil
from retrying import retry

from data_refinery_common.job_lookup import (
    ProcessorPipeline,
    determine_processor_pipeline,
    determine_ram_amount,
)
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
from data_refinery_common.utils import get_env_variable, get_instance_id

logger = get_and_configure_logger(__name__)
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")

CURRENT_JOB = None


def signal_handler(sig, frame):
    """Signal Handler, works for both SIGTERM and SIGINT"""
    global CURRENT_JOB
    if CURRENT_JOB:
        CURRENT_JOB.success = False
        CURRENT_JOB.end_time = timezone.now()
        CURRENT_JOB.num_retries = CURRENT_JOB.num_retries - 1
        CURRENT_JOB.failure_reason = "Interruped by SIGTERM/SIGINT: " + str(sig)
        CURRENT_JOB.save()

    sys.exit(0)


def start_job(job_id: int) -> DownloaderJob:
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

    # This job should not have been started.
    if job.start_time is not None:
        logger.error("This downloader job has already been started!!!", downloader_job=job.id)
        raise Exception("downloaders.start_job called on a job that has already been started!")

    # Set up the SIGTERM handler so we can appropriately handle being interrupted.
    # (`docker stop` uses SIGTERM, not SIGINT.)
    # (however, Nomad sends an SIGINT so catch both.)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    job.worker_id = worker_id
    job.worker_version = SYSTEM_VERSION
    job.start_time = timezone.now()
    job.save()

    needs_downloading = any(
        original_file.needs_downloading() for original_file in job.original_files.all()
    )

    if not needs_downloading:
        logger.error(
            ("No files associated with this job need to be downloaded! Aborting!"), job_id=job.id
        )
        job.start_time = timezone.now()
        job.failure_reason = "Was told to redownload file(s) that are already downloaded!"
        job.success = False
        job.no_retry = True
        job.end_time = timezone.now()
        job.save()
        sys.exit(0)

    global CURRENT_JOB
    CURRENT_JOB = job

    return job


def end_downloader_job(job: DownloaderJob, success: bool):
    """
    Record in the database that this job has completed.
    """
    if success:
        logger.debug("Downloader Job completed successfully.", downloader_job=job.id)
    else:
        # Should be set by now, but make sure.
        success = False
        file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
        for file_assoc in file_assocs:
            file_assoc.original_file.delete_local_file()
            file_assoc.original_file.is_downloaded = False
            file_assoc.original_file.save()

        if not job.failure_reason:
            logger.error(
                "Downloader job failed without having failure_reason set. FIX ME!!!!!!!!",
                downloader_job=job.id,
                downloader_task=job.downloader_task,
            )
        else:
            logger.info(
                "Downloader job failed!",
                downloader_job=job.id,
                downloader_task=job.downloader_task,
                failure_reason=job.failure_reason,
            )

    job.success = success
    job.end_time = timezone.now()
    job.save()

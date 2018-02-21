from django.utils import timezone
from django.db import transaction
from retrying import retry
from typing import List, Dict

from data_refinery_common.utils import get_worker_id
from data_refinery_common.models import (
    Batch,
    BatchStatuses,
    DownloaderJob,
    ProcessorJob,
    Experiment,
    Sample,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_workers._version import __version__


logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "downloader_job_"


def start_job(job_id: int) -> DownloaderJob:
    """Record in the database that this job is being started.

    Retrieves the job from the database and returns it after marking
    it as started.
    """
    logger.info("Starting Downloader Job.", downloader_job=job_id)
    try:
        job = DownloaderJob.objects.get(id=job_id)
    except DownloaderJob.DoesNotExist:
        logger.error("Cannot find downloader job record.", downloader_job=job_id)
        raise

    job.worker_id = get_worker_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    return job

def end_downloader_job(job: DownloaderJob, success: bool):
    """
    Record in the database that this job has completed.
    """

    job.success = success
    job.end_time = timezone.now()
    job.save()

def create_processor_jobs_for_original_files(original_files: List[OriginalFile]):
    """
    Create a processor jobs queue a processor task for samples related to an experiment.
    """

    # Iterate over all of our samples.
    # If we have raw, send it to the correct processor.
    # Else, treat it as a "NO-OP"
    for original_file in original_files:

        processor_job = ProcessorJob()

        if not original_file.has_raw:
            processor_job.pipeline_applied = ProcessorPipeline.NO_OP
        else:
            processor_job.pipeline_applied = "AFFY_TO_PCL"

        # Save the Job and create the association
        processor_job.save()
        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job
        assoc.save()

        send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)

def create_processor_job_for_original_files(original_files: List[OriginalFile]):
    """
    Create a processor job and queue a processor task for sample related to an experiment.

    """
    original_file = original_files[0]

    # This is a paired read. Make sure the other one is downloaded, this start the job
    if '_' in original_file.file_name:
        split = original_file.file_name.split('_')
        if '1' in split[1]:
            other_file = OriginalFile.objects.get(source_filename='_'.join([split[0], split[1].replace('1', '2')]))
        else:
            other_file = OriginalFile.objects.get(source_filename='_'.join([split[0], split[1].replace('2', '1')]))
        if not other_file.is_downloaded:
            logger.info("Need other file to download before starting paired read Salmon.")
            return
        else:
            processor_job = ProcessorJob()
            processor_job.pipeline_applied = "SALMON"
            processor_job.save()

            assoc1 = ProcessorJobOriginalFileAssociation()
            assoc1.original_file = original_file
            assoc1.processor_job = processor_job
            assoc1.save()

            assoc2 = ProcessorJobOriginalFileAssociation()
            assoc2.original_file = other_file
            assoc2.processor_job = processor_job
            assoc2.save()

    # This is a single read. Let's rock now.
    else:
        # Only one file to download, start the job now.
        processor_job.save()
        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job
        assoc.save()

    send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)


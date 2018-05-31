from django.utils import timezone
from django.db import transaction
from retrying import retry
from typing import List, Dict

from data_refinery_common.utils import get_worker_id
from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    Experiment,
    Sample,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_common.job_lookup import ProcessorPipeline, determine_processor_pipeline
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
    for original_file in original_files:
        # sample_object = Sample.objects.filter(original_file=original_file).first()
        # Might work?
        sample_object = original_file.samples.first().accession_code

        processor_job = ProcessorJob()
        processor_job.pipeline_applied = determine_processor_pipeline(sample_object)
        processor_job.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job
        assoc.save()

        send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)


def create_processor_job_for_original_files(original_files: List[OriginalFile],
                                            sample_object: Sample=None):
    """
    Create a processor job and queue a processor task for sample related to an experiment.

    """
    if not sample_object:
        # XXX: do this right, but for now I wanna keep moving
        # I should probably get the samples for each one and make sure they are the same sample.
        # Also consider what happens if there isn't one? That's a pretty BFD
        # Actual comment:
        # For anything that has raw data there should only be one Sample per OriginalFile
        # sample_object = Sample.objects.filter(original_file=original_files[0]).first()
        # Might work?
        sample_object = original_files[0].samples.first().accession_code

    processor_job = ProcessorJob()
    processor_job.pipeline_applied = determine_processor_pipeline(sample_object)
    processor_job.save()
    for original_file in original_files:
        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job
        assoc.save()

    send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)

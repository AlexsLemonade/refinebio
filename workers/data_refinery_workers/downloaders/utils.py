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
from data_refinery_common.job_lookup import ProcessorPipeline, determine_processor_pipeline, determine_ram_amount
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
    if success:
        logger.info("Downloader Job completed successfully.",
                    downloader_job=job.id)
    else:
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


def create_processor_jobs_for_original_files(original_files: List[OriginalFile],
                                             downloader_job: DownloaderJob=None):
    """
    Create a processor jobs and queue a processor task for samples related to an experiment.
    """
    for original_file in original_files:
        sample_object = original_file.samples.first()

        if original_file.filename[-4:] == ".xml":
            logger.info("Skipping useless file",
                file=original_file.id,
                filename=original_file.filename)
            continue

        # We NO_OP processed data. It's what we do.
        if '.processed' in original_file.source_url:
            pipeline_to_apply = ProcessorPipeline.NO_OP
        else:
            pipeline_to_apply = determine_processor_pipeline(sample_object, original_file)

        if pipeline_to_apply == ProcessorPipeline.NONE:
            logger.info("No valid processor pipeline found to apply to sample.",
                        sample=sample_object.id,
                        original_file=original_files[0].id)
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
                logger.info("Queuing processor job.",
                            processor_job=processor_job.id,
                            original_file=original_file.id,
                            downloader_job=downloader_job.id)
            else:
                logger.info("Queuing processor job.",
                            processor_job=processor_job.id,
                            original_file=original_file.id)

            send_job(pipeline_to_apply, processor_job)


def create_processor_job_for_original_files(original_files: List[OriginalFile],
                                            downloader_job: DownloaderJob=None):
    """
    Create a processor job and queue a processor task for sample related to an experiment.

    """

    # For anything that has raw data there should only be one Sample per OriginalFile
    sample_object = original_files[0].samples.first()
    pipeline_to_apply = determine_processor_pipeline(sample_object, original_files[0])

    if pipeline_to_apply == ProcessorPipeline.NONE:
        logger.info("No valid processor pipeline found to apply to sample.",
                    sample=sample_object.id,
                    original_file=original_files[0].id)
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

        if downloader_job:
            logger.info("Queuing processor job.",
                        processor_job=processor_job.id,
                        downloader_job=downloader_job.id)
        else:
            logger.info("Queuing processor job.",
                        processor_job=processor_job.id)

        send_job(pipeline_to_apply, processor_job)

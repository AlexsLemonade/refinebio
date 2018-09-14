from django.db import transaction
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
from data_refinery_common.utils import get_instance_id


logger = get_and_configure_logger(__name__)
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")
# TODO: extend this list.
BLACKLISTED_EXTENSIONS = ["xml", "chp", "exp"]


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

    # This job should not have been started.
    if job.start_time is not None:
        logger.error("This downloader job has already been started!!!", downloader_job=job.id)
        raise Exception("downloaders.start_job called on a job that has already been started!")

    job.worker_id = get_instance_id()
    job.worker_version = SYSTEM_VERSION
    job.start_time = timezone.now()
    job.save()

    return job


def end_downloader_job(job: DownloaderJob, success: bool):
    """
    Record in the database that this job has completed.
    """
    if success:
        logger.debug("Downloader Job completed successfully.",
                    downloader_job=job.id)
    else:
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

    return OriginalFile


def create_processor_jobs_for_original_files(original_files: List[OriginalFile],
                                             downloader_job: DownloaderJob=None):
    """
    Create a processor jobs and queue a processor task for samples related to an experiment.
    """
    for original_file in original_files:
        sample_object = original_file.samples.first()

        if not delete_if_blacklisted(original_file):
            continue

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

        send_job(pipeline_to_apply, processor_job)

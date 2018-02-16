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
)
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.models.new_models import Experiment, Sample, ExperimentAnnotation, ExperimentSampleAssociation, OriginalFile, ProcessorJobOriginalFileAssociation
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


# def end_job(job: DownloaderJob, batches: Batch, success: bool):
#     """Record in the database that this job has completed.

#     Create a processor job and queue a processor task for each batch
#     if the job was successful.
#     """
#     @retry(stop_max_attempt_number=3)
#     def save_batch_create_job(batch):
#         batch.status = BatchStatuses.DOWNLOADED.value
#         batch.save()

#         # TEMPORARY for Jackie's grant:
#         if batch.pipeline_required != ProcessorPipeline.NONE.value:
#             logger.debug("Creating processor job for Batch.",
#                          downloader_job=job.id,
#                          batch=batch.id)
#             with transaction.atomic():
#                 processor_job = ProcessorJob.create_job_and_relationships(
#                     batches=[batch], pipeline_applied=batch.pipeline_required)
#             return processor_job
#         else:
#             logger.debug("Not queuing a processor job for batch.",
#                          downloader_job=job.id,
#                          batch=batch.id)
#             return None

#     @retry(stop_max_attempt_number=3)
#     def queue_task(processor_job, batch):
#         if batch.pipeline_required in ProcessorPipeline.__members__:
#             send_job(ProcessorPipeline[batch.pipeline_required], processor_job.id)
#             logger.info("Queuing processor job.",
#                         downloader_job=job.id,
#                         processor_job=processor_job.id,
#                         batch=batch.id)
#             return True
#         else:
#             failure_template = "Could not find Processor Pipeline {} in the lookup."
#             failure_message = failure_template.format(batch.pipeline_required)
#             logger.error(failure_message, downloader_job=job.id, batch=batch.id)
#             processor_job.failure_reason = failure_message
#             processor_job.success = False
#             processor_job.retried = True
#             processor_job.save()
#             return False

#     if success:
#         for batch in batches:
#             processor_job = save_batch_create_job(batch)
#             if batch.pipeline_required != ProcessorPipeline.NONE.value:
#                 try:
#                     success = queue_task(processor_job, batch)
#                 except:
#                     logger.exception("Could not queue processor job task.")
#                     # If the task doesn't get sent we don't want the
#                     # processor_job to be left floating
#                     processor_job.delete()

#                     success = False
#                     job.failure_message = "Could not queue processor job task."

#                 if success:
#                     logger.info("Downloader job completed successfully.", downloader_job=job.id)

#     # Check to make sure job didn't end because of missing batches or files.
#     if len(batches) > 0 and len(batches[0].files) > 0:
#         # Clean up temp directory to free up local disk space.
#         batches[0].files[0].remove_temp_directory(JOB_DIR_PREFIX + str(job.id))

#     job.success = success
#     job.end_time = timezone.now()
#     job.save()

def end_downloader_job(job: DownloaderJob, success: bool):
    """
    Record in the database that this job has completed.
    """

    job.success = success
    job.end_time = timezone.now()
    job.save()

# def create_processor_jobs_for_experiment(experiment: Experiment):
#     """
#     Create a processor job and queue a processor task for sample related to an experiment.
#     """

#     relations = ExperimentSampleAssociation.objects.filter(experiment=experiment)
#     samples = Sample.objects.filter(id__in=relations.values('sample_id'))

#     # Iterate over all of our samples.
#     # If we have raw, send it to the correct processor.
#     # Else, treat it as a "NO-OP"
#     for sample in samples:

#         processor_job = ProcessorJob()
#         if not sample.has_raw:
#             processor_job.pipeline_applied = "NO_OP"
#         else:
#             if experiment.source_database == "ARRAY_EXPRESS":
#                 processor_job.pipeline_applied = "AFFY_TO_PCL"
#             # TODO: Salmon!

#         # Save the Job and create the association
#         processor_job.save()
#         assoc = ProcessorJobSampleAssociation()
#         assoc.sample = sample
#         assoc.processor_job = processor_job
#         assoc.save()

#         # Send the job to Nomad
#         send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)

def create_processor_jobs_for_original_files(original_files: List[OriginalFile]):
    """
    Create a processor job and queue a processor task for sample related to an experiment.
    
    TODO: Refactor to support multifile Salmon?

    """

    processor_job = ProcessorJob()

    # Iterate over all of our samples.
    # If we have raw, send it to the correct processor.
    # Else, treat it as a "NO-OP"
    for original_file in original_files:

        if original_file.sample:
            experiment = ExperimentSampleAssociation.objects.filter(sample=original_file.sample)[0].experiment

        if experiment.source_database == "ARRAY_EXPRESS":
            if not original_file.has_raw:
                processor_job.pipeline_applied = ProcessorPipeline.NO_OP
            else:
                processor_job.pipeline_applied = "AFFY_TO_PCL"
        elif experiment.source_database == "SRA":
            processor_job.pipeline_applied = "SALMON" 

        # Save the Job and create the association
        processor_job.save()
        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job
        assoc.save()

    # Send the job to Nomad
    send_job(ProcessorPipeline[processor_job.pipeline_applied], processor_job.id)


    # @retry(stop_max_attempt_number=3)
    # def save_batch_create_job(batch):
    #     batch.status = BatchStatuses.DOWNLOADED.value
    #     batch.save()

    #     # TEMPORARY for Jackie's grant:
    #     if batch.pipeline_required != ProcessorPipeline.NONE.value:
    #         logger.debug("Creating processor job for Batch.",
    #                      downloader_job=job.id,
    #                      batch=batch.id)
    #         with transaction.atomic():
    #             processor_job = ProcessorJob.create_job_and_relationships(
    #                 batches=[batch], pipeline_applied=batch.pipeline_required)
    #         return processor_job
    #     else:
    #         logger.debug("Not queuing a processor job for batch.",
    #                      downloader_job=job.id,
    #                      batch=batch.id)
    #         return None

    # @retry(stop_max_attempt_number=3)
    # def queue_task(processor_job, batch):
    #     if batch.pipeline_required in ProcessorPipeline.__members__:
    #         send_job(ProcessorPipeline[batch.pipeline_required], processor_job.id)
    #         logger.info("Queuing processor job.",
    #                     downloader_job=job.id,
    #                     processor_job=processor_job.id,
    #                     batch=batch.id)
    #         return True
    #     else:
    #         failure_template = "Could not find Processor Pipeline {} in the lookup."
    #         failure_message = failure_template.format(batch.pipeline_required)
    #         logger.error(failure_message, downloader_job=job.id, batch=batch.id)
    #         processor_job.failure_reason = failure_message
    #         processor_job.success = False
    #         processor_job.retried = True
    #         processor_job.save()
    #         return False

    # for batch in batches:
    #     processor_job = save_batch_create_job(batch)
    #     if batch.pipeline_required != ProcessorPipeline.NONE.value:
    #         try:
    #             success = queue_task(processor_job, batch)
    #         except:
    #             logger.exception("Could not queue processor job task.")
    #             # If the task doesn't get sent we don't want the
    #             # processor_job to be left floating
    #             processor_job.delete()

    #             success = False
    #             job.failure_message = "Could not queue processor job task."

    #         if success:
    #             logger.info("Downloader job completed successfully.", downloader_job=job.id)

    # # Check to make sure job didn't end because of missing batches or files.
    # if len(batches) > 0 and len(batches[0].files) > 0:
    #     # Clean up temp directory to free up local disk space.
    #     batches[0].files[0].remove_temp_directory(JOB_DIR_PREFIX + str(job.id))


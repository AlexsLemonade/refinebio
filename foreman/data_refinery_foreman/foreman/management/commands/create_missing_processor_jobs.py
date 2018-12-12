"""
This command will run the Foreman's main function monitor_jobs.
This will cause the Foreman to check for a number of different
failures for both the DownloaderJobs and ProcessorJobs and requeue
those jobs it detects as failed.
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    OriginalFile,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_common.utils import is_file_rnaseq
from data_refinery_foreman.foreman.shepherd import create_missing_processor_jobs


def create_processor_jobs_for_original_files(original_files: List[OriginalFile],
                                             volume_index: int):
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
            processor_job.volume_index = volume_index
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
                                            volume_index: int):
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
        processor_job.volume_index = volume_index
        processor_job.save()
        for original_file in original_files:
            assoc = ProcessorJobOriginalFileAssociation()
            assoc.original_file = original_file
            assoc.processor_job = processor_job
            assoc.save()
        logger.debug("Queuing processor job.",
                     processor_job=processor_job.id)
        send_job(pipeline_to_apply, processor_job)


def find_volume_index_for_dl_job(job: DownloaderJob) -> int:
    pjs = ProcessorJob.objects.filter(worker_id=job.worker_id)

    if pjs:
        retrun pjs.first().volume_index
    else:
        logger.warn(
            "Could not determine what volume index a downloader job was run on",
            downloader_job=job.id
        )
        raise ValueError()


class Command(BaseCommand):
    def handle(self, *args, **options):
        for dl_job in DownloaderJob.objects.filter(success='t').all():
            original_files = dl_job.original_files

            if original_files.count() > 0:
                if is_file_rnaseq(original_files.first().filename):
                    try:
                        find_volume_index_for_dl_job(dl_job)
                        create_processor_job_for_original_files(original_files, volume_index)
                    except:
                        # Already logged.
                        pass
                else:
                    non_archive_files = []
                    for og_file in dl_job.original_files:
                        if not og_file.is_archive:
                            non_archive_files.append(og_file)

                    if non_archive_files:
                        try:
                            find_volume_index_for_dl_job(dl_job)
                            create_processor_jobs_for_original_files(
                                non_archive_files,
                                volume_index
                            )
                        except:
                            # Already logged.
                            pass

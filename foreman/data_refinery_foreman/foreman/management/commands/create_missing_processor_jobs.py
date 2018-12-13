"""This command will find successful downloader jobs who never
spawned a processor job and will create processor jobs for them. This
is most likely to be a temporary patch until we can figure out why
downloader jobs are failing to create processor jobs.
"""

from typing import List

from django.core.management.base import BaseCommand
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    OriginalFile,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_common.job_lookup import is_file_rnaseq


logger = get_and_configure_logger(__name__)


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
        return pjs.first().volume_index
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
                first_file=original_files.first()
                if is_file_rnaseq(first_file.filename):
                    try:
                        sample_object = first_file.samples.first()
                        logger.info(("Found a downloader job that didn't make a processor"
                                     " job for an RNA-Seq sample."),
                                    downloader_job=dl_job.id,
                                    file_name=first_file.filename,
                                    original_file_id=first_file.id,
                                    sample_accession=sample_object.accession_code,
                                    sample_id=sample_object.id
                        )
                        find_volume_index_for_dl_job(dl_job)
                        create_processor_job_for_original_files(original_files, volume_index)
                    except:
                        # Already logged.
                        pass
                else:
                    # When the GEO downloader extracts files it
                    # downloads, it doesn't link them back to the
                    # downloader job, but they do get linked to the
                    # sample.
                    all_original_files = first_file.samples.first().original_files.all()
                    non_archive_files = []
                    for og_file in all_original_files:
                        if not og_file.is_archive:
                            non_archive_files.append(og_file)

                    if non_archive_files:
                        try:
                            sample_object = non_archive_files[0].samples.first()
                            logger.info(("Found a downloader job that didn't make a processor"
                                         " job for an Microarray sample."),
                                        downloader_job=dl_job.id,
                                        num_files_in_dl_job=len(non_archive_files),
                                        sample_accession=sample_object.accession_code,
                                        sample_id=sample_object.id
                            )

                            find_volume_index_for_dl_job(dl_job)
                            create_processor_jobs_for_original_files(
                                non_archive_files,
                                volume_index
                            )
                        except:
                            # Already logged.
                            pass

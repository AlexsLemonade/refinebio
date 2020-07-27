"""This command will find successful downloader jobs who never
spawned a processor job and will create processor jobs for them. This
is most likely to be a temporary patch until we can figure out why
downloader jobs are failing to create processor jobs.
"""


from django.core.management.base import BaseCommand

from data_refinery_common.job_lookup import (
    ProcessorPipeline,
    determine_processor_pipeline,
    determine_ram_amount,
    is_file_rnaseq,
)
from data_refinery_common.job_management import (
    create_processor_job_for_original_files,
    create_processor_jobs_for_original_files,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
)

logger = get_and_configure_logger(__name__)


def find_volume_index_for_dl_job(job: DownloaderJob) -> int:
    pjs = ProcessorJob.objects.filter(worker_id=job.worker_id)

    if pjs:
        return pjs.first().volume_index
    else:
        logger.warn(
            "Could not determine what volume index a downloader job was run on",
            downloader_job=job.id,
        )
        raise ValueError()


class Command(BaseCommand):
    def handle(self, *args, **options):
        for dl_job in DownloaderJob.objects.filter(success="t").all():
            original_files = dl_job.original_files.all()

            if original_files.count() > 0:
                first_file = original_files.first()
                if (
                    is_file_rnaseq(first_file.filename)
                    and first_file.processor_jobs.all().count() == 0
                ):
                    try:
                        sample_object = first_file.samples.first()
                        logger.info(
                            (
                                "Found a downloader job that didn't make a processor"
                                " job for an RNA-Seq sample."
                            ),
                            downloader_job=dl_job.id,
                            file_name=first_file.filename,
                            original_file_id=first_file.id,
                            sample_accession=sample_object.accession_code,
                            sample_id=sample_object.id,
                        )
                        volume_index = find_volume_index_for_dl_job(dl_job)
                        create_processor_job_for_original_files(
                            original_files, dl_job, volume_index
                        )
                    except Exception:
                        # Already logged.
                        pass
                else:
                    # When the GEO downloader extracts files it
                    # downloads, it doesn't link them back to the
                    # downloader job, but they do get linked to the
                    # sample.
                    samples_in_job = set()
                    for og_file in original_files:
                        for sample in og_file.samples.all():
                            samples_in_job.add(sample)

                    for sample in samples_in_job:
                        all_original_files = sample.original_files.all()
                        files_for_sample = []
                        for og_file in all_original_files:
                            if not og_file.is_archive and og_file.processor_jobs.all().count() == 0:
                                files_for_sample.append(og_file)

                        if files_for_sample:
                            try:
                                logger.info(
                                    (
                                        "Found a downloader job that didn't make a processor"
                                        " job for an Microarray sample."
                                    ),
                                    downloader_job=dl_job.id,
                                    sample_accession=sample.accession_code,
                                    sample_id=sample.id,
                                )

                                volume_index = find_volume_index_for_dl_job(dl_job)
                                create_processor_jobs_for_original_files(
                                    files_for_sample, dl_job, volume_index
                                )
                            except Exception:
                                # Already logged.
                                pass

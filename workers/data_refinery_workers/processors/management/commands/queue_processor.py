"""This command is intended for development purposes.
It creates the database records necessary for a processor job to
run and queues one. It assumes that the file
/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072_CD_colon_active_2.CEL
exists.
The easiest way to run this is with the tester.sh script.
(Changing queue_downloader to queue_processor.)
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob
)
from data_refinery_workers.task_runner import app
from data_refinery_common.job_lookup import PROCESSOR_PIPELINE_LOOKUP
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        # Create all the dummy data that would have been created
        # before a processor job could have been generated.
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="SRA",
            pipeline_required="SALMON",
            platform_accession_code="IlluminaHiSeq2500",
            experiment_accession_code="DRX014494",
            experiment_title="It doesn't really matter.",
            organism_id=10090,
            organism_name="ARABIDOPSIS THALIANA",
            release_date="2015-05-03",
            last_uploaded_date="2015-06-19",
            status=BatchStatuses.DOWNLOADED.value,
        )
        batch.save()

        file1 = File(
            size_in_bytes=967794560,
            raw_format="fastq.gz",
            processed_format="tar.gz",
            name="DRR016125_1.fastq.gz",
            internal_location="IlluminaHiSeq2500/SALMON",
            download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_1.fastq.gz",  # noqa
            batch=batch
        )
        file1.save()

        file2 = File(
            size_in_bytes=1001146319,
            raw_format="fastq.gz",
            processed_format="tar.gz",
            name="DRR016125_2.fastq.gz",
            internal_location="IlluminaHiSeq2500/SALMON",
            download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_2.fastq.gz",  # noqa
            batch=batch
        )
        file2.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        processor_task = PROCESSOR_PIPELINE_LOOKUP[batch.pipeline_required]
        app.send_task(processor_task, args=[processor_job.id])
        logger.info("Processor Job queued.", processor_job=processor_job.id)

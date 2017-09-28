"""This command is intended for development purposes.
It creates the database records necessary for a processor job to
run and queues one. It assumes that the file
/home/user/data_store/raw/A-AFFY-141/AFFY_TO_PCL/GSM1426072_CD_colon_active_2.CEL
exists.
The easiest way to run this is with the tester.sh script.
(Changing queue_downloader to queue_processor.)
"""

from django.core.management.base import BaseCommand
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    ProcessorJob
)
from data_refinery_workers.task_runner import app
from data_refinery_common.job_lookup import PROCESSOR_PIPELINE_LOOKUP


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        # Create all the dummy data that would have been created
        # before a processor job could have been generated.
        survey_job = SurveyJob(
            source_type="SRA"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="SRA",
            size_in_bytes=2214725074,
            raw_format="fastq",
            processed_format="sf",
            pipeline_required="SALMON",
            platform_accession_code="IlluminaHiSeq2500",
            experiment_accession_code="PRJEB5018",
            experiment_title="It doesn't really matter.",
            name="ERR1680082_1.fastq",
            internal_location="IlluminaHiSeq2500/SALMON",
            organism_id=10090,
            organism_name="MUS MUSCULUS",
            release_date="2014-03-25",
            last_uploaded_date="2016-05-20",
            status=BatchStatuses.NEW.value,
            download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/ERR1680082/ERR1680082_1.fastq.gz"  # noqa
        )
        batch.save()

        batch2 = Batch(
            survey_job=survey_job,
            source_type="SRA",
            size_in_bytes=2214725074,
            download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/ERR1680082/ERR1680082_2.fastq.gz",  # noqa
            raw_format="fastq",
            processed_format="sf",
            pipeline_required="SALMON",
            platform_accession_code="IlluminaHiSeq2500",
            experiment_accession_code="PRJEB5018",
            experiment_title="It doesn't really matter.",
            name="ERR1680082_2.fastq",
            internal_location="IlluminaHiSeq2500/SALMON",
            organism_id=10090,
            organism_name="MUS MUSCULUS",
            release_date="2014-03-25",
            last_uploaded_date="2016-05-20",
            status=BatchStatuses.NEW.value
        )
        batch2.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch, batch2])
        logger.info("Queuing a processor job.")
        processor_task = PROCESSOR_PIPELINE_LOOKUP[batch.pipeline_required]
        app.send_task(processor_task, args=[processor_job.id])

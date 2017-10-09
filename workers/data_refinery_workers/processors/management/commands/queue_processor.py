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
        survey_job = SurveyJob(
            source_type="ARRAY_EXPRESS"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="ARRAY_EXPRESS",
            size_in_bytes=0,
            download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",  # noqa
            raw_format="CEL",
            processed_format="PCL",
            pipeline_required="NO_OP",
            platform_accession_code="A-AFFY-141",
            experiment_accession_code="E-GEOD-59071",
            experiment_title="It doesn't really matter.",
            name="GSM1426074_CD_colon_active_4.CEL",
            internal_location="A-AFFY-141/AFFY_TO_PCL",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        processor_task = PROCESSOR_PIPELINE_LOOKUP[batch.pipeline_required]
        app.send_task(processor_task, args=[processor_job.id])
        logger.info("Processor Job queued.", processor_job=processor_job.id)

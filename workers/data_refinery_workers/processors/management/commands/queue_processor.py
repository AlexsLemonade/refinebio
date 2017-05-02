"""This command is intended for development purposes.
It creates the database records necessary for a processor job to
run and queues one. It assumes that the file
/home/user/data_store/raw/A-AFFY-1/MICRO_ARRAY_TO_PCL/E-MTAB-3050.raw.1.zip
exists.
The easiest way to run this is with the tester.sh script.
(Changing queue_downloader to queue_processor.)"""

from django.core.management.base import BaseCommand
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    ProcessorJob
)
from data_refinery_workers.processors.array_express \
    import process_array_express


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
            download_url="http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3050/E-MTAB-3050.raw.1.zip",  # noqa
            raw_format="MICRO_ARRAY",
            processed_format="PCL",
            pipeline_required="MICRO_ARRAY_TO_PCL",
            accession_code="A-AFFY-1",
            internal_location="A-AFFY-1/MICRO_ARRAY_TO_PCL/",
            organism=1,
            status=BatchStatuses.NEW.value
        )
        batch.save()

        processor_job = ProcessorJob(batch=batch)
        processor_job.save()
        logger.info("Queuing a processor job.")
        process_array_express.delay(processor_job.id)

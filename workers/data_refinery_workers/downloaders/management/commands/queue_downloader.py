"""This command is intended for development purposes.
It creates the database records necessary for a downloader job to
run and queues one.
The easiest way to run this is with the tester.sh script."""

from django.core.management.base import BaseCommand
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    DownloaderJob,
    DownloaderJobsToBatches
)
from data_refinery_workers.downloaders.array_express \
    import download_array_express


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        # Create all the dummy data that would have been created
        # before a downloader job could have been generated.
        survey_job = SurveyJob(
            source_type="ARRAY_EXPRESS"
        )
        survey_job.save()

        batch = Batch(
            survey_job=survey_job,
            source_type="ARRAY_EXPRESS",
            size_in_bytes=0,
            download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip/GSM1426072_CD_colon_active_2.CEL",  # noqa
            raw_format="CEL",
            processed_format="PCL",
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-1",
            experiment_accession_code="E-MTAB-3050",
            experiment_title="It doesn't really matter.",
            name="CE12314",
            internal_location="A-AFFY-1/MICRO_ARRAY_TO_PCL/",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()

        downloader_job = DownloaderJob()
        downloader_job.save()
        downloader_job_to_batch = DownloaderJobsToBatches(batch=batch,
                                                          DownloaderJob=downloader_job)
        downloader_job_to_batch.save()
        logger.info("Queuing a task.")
        download_array_express.delay(downloader_job.id)

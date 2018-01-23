"""This command is intended for development purposes.
It creates the database records necessary for a downloader job to
run and queues one.
The easiest way to run this is with the tester.sh script."""

from django.core.management.base import BaseCommand
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    DownloaderJob
)


logger = get_and_configure_logger(__name__)


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
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-141",
            experiment_accession_code="E-GEOD-59071",
            experiment_title="It doesn't really matter.",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()

        file = File(
            batch=batch,
            size_in_bytes=0,
            download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",  # noqa
            raw_format="CEL",
            processed_format="PCL",
            name="GSM1426072_CD_colon_active_2.CEL",
            internal_location="A-AFFY-141/AFFY_TO_PCL"
        )
        file.save()

        downloader_job = DownloaderJob.create_job_and_relationships(batches=[batch])
        send_job(Downloaders["ARRAY_EXPRESS"], downloader_job.id)

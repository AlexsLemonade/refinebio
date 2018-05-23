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
    Organism,
    Experiment,
    Sample,
    OriginalFile,
    ExperimentSampleAssociation,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation
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

        organism = Organism.get_object_for_name("HOMO_SAPIENS")

        # Remove already existing objects so these can be re-created.
        Experiment.objects.filter(accession_code="E-GEOD-59071").delete()
        Sample.objects.filter(accession_code="GSM1426072").delete()

        experiment_object = Experiment()
        experiment_object.accession_code = "E-GEOD-59071"
        experiment_object.source_url = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-59071/"
        experiment_object.source_database = "ARRAY_EXPRESS"
        experiment_object.name = "It doesn't really matter."
        experiment_object.description = "It doesn't really matter."
        experiment_object.platform_name = "Some platform name..."
        experiment_object.platform_accession_code = "A-AFFY-141"
        experiment_object.source_first_published = "2017-05-05"
        experiment_object.source_last_modified = "2017-05-05"
        experiment_object.save()

        sample_object = Sample()
        sample_object.accession_code = "GSM1426072"
        sample_object.source_archive_url = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-59071/samples"
        sample_object.organism = organism
        sample_object.save()

        original_file = OriginalFile()
        original_file.sample = sample_object
        original_file.source_filename = "GSM1426072_CD_colon_active_2.CEL"
        original_file.source_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        original_file.is_downloaded = False
        original_file.is_archive = True
        original_file.has_raw = True
        original_file.save()

        association = ExperimentSampleAssociation()
        association.experiment = experiment_object
        association.sample = sample_object
        association.save()

        downloader_job = DownloaderJob()
        downloader_job.downloader_task = Downloaders.ARRAY_EXPRESS.value
        downloader_job.accession_code = experiment_object.accession_code
        downloader_job.save()

        asoc = DownloaderJobOriginalFileAssociation()
        asoc.downloader_job = downloader_job
        asoc.original_file = original_file
        asoc.save()

        logger.info("Queuing downloader job.",
                    survey_job=survey_job.id,
                    downloader_job=downloader_job.id)
        send_job(Downloaders[downloader_job.downloader_task], downloader_job.id)

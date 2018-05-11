import time
import json
from django.test import TransactionTestCase, tag
from datetime import timedelta, datetime
from django.utils import timezone
from unittest.mock import patch, Mock
from data_refinery_common.models import (
    Organism,
    SurveyJob,
    SurveyJobKeyValue,
    DownloaderJob,
    ProcessorJob,
    Sample
)
from data_refinery_foreman.surveyor import surveyor

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
module_logger = logging.getLogger("data_refinery_foreman.surveyor.array_express")


LOOP_TIME = 5  # seconds
MAX_WAIT_TIME = timedelta(minutes=5)

def wait_for_job(job, job_class: type, start_time: datetime):
    """Monitors the `job_class` table for when `job` is done."""
    job = job_class.objects.filter(id=job.id).get()
    while job.success is None and timezone.now() - start_time < MAX_WAIT_TIME:
        logger.info("Still polling the %s.",
                    job_class.__name__)
        time.sleep(LOOP_TIME)
        job = job_class.objects.filter(id=job.id).get()

    return job


def mock_get_sample(accession_code: str):
    """Mock out Sample.objects.get to prevent too much work from being done.

    There are a lot of Samples on this experiment, but we only want to
    do one from each zip file. To accomplish this we make every other
    Sample in the experiment appear to already exist so it is not
    created and processed.
    """
    if accession_code == "GSM1109016" or accession_code == "GSM1108516":
        raise Sample.DoesNotExist

    # This return value isn't actually used so it doesn't matter what
    # it is. The call to this function is just checking for the
    # existence of the Sample.
    return None

def mock_logger_info(message: str, *args, **kwargs):
    """Silence the log messages we're forcing on purpose.

    Since we fake out having a ton of samples in mock_get_sample, the
    surveyor logs a bunch of messages about them already
    existing. However we're expecting this and it blows up the logs
    for the test. Therefore if the message is about that, don't
    log.
    """
    if "skipping object creation" in message:
        return None

    module_logger.info(message, *args, **kwargs)


# TransactionTestCase makes database calls complete before the test
# ends.  Otherwise the workers wouldn't actually be able to find the
# job in the database cause it'd be stuck in a transaction.
class NoOpEndToEndTestCase(TransactionTestCase):
    @tag("slow")
    @patch('data_refinery_foreman.surveyor.array_express.logger')
    @patch('data_refinery_common.models.Sample.objects.get')
    def test_no_op(self, mocked_get_query, mocked_logger):
        """Survey, download, then process an experiment we know is NO_OP."""
        mocked_get_query.side_effect = mock_get_sample
        mocked_logger.side_effect = mock_logger_info

        # Prevent a call being made to NCBI's API to determine
        # organism name/id.
        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        survey_job = surveyor.survey_ae_experiment("E-GEOD-45547")

        self.assertTrue(survey_job.success)

        downloader_jobs = DownloaderJob.objects.all()
        logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
        start_time = timezone.now()
        for downloader_job in downloader_jobs:
            downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time)
            self.assertTrue(downloader_job.success)

        processor_jobs = ProcessorJob.objects.all()
        logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
        start_time = timezone.now()
        for processor_job in processor_jobs:
            processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
            self.assertTrue(processor_job.success)

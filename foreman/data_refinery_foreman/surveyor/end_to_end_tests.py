import time
from django.test import TransactionTestCase
from data_refinery_models.models import (
    SurveyJob,
    SurveyJobKeyValue,
    Batch,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_foreman.surveyor import surveyor

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


LOOP_TIME = 5  # seconds


def wait_for_job(job, job_class: type):
    """Monitors the `job_class` table for when `job` is done."""
    job = job_class.objects.filter(id=job.id).get()
    while job.success is None:
        logger.info("Still polling the %s.",
                    job_class.__name__)
        time.sleep(LOOP_TIME)
        job = job_class.objects.filter(id=job.id).get()

    return job


class ScanUpcEndToEndTestCase(TransactionTestCase):
    def test_calls_survey(self):
        """If source_type is supported calls the appropriate survey method."""
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value="E-GEOD-22166")
        key_value_pair.save()

        surveyor.run_job(survey_job)
        logger.info("Started Survey Job %s, waiting for it to complete.", survey_job.id)
        survey_job = wait_for_job(survey_job, SurveyJob)
        self.assertTrue(survey_job.success)

        batch = Batch.objects.all()[0]
        batch = Batch.objects.filter(survey_job=survey_job).get()

        downloader_job = batch.downloaderjob_set.get()
        logger.info("Survey Job finished, waiting for Downloader Job %d to complete.",
                    downloader_job.id)
        downloader_job = wait_for_job(downloader_job, DownloaderJob)
        self.assertTrue(downloader_job.success)

        processor_job = batch.processorjob_set.get()
        logger.info("Downloader Job finished, waiting for processor Job %d to complete.",
                    processor_job.id)
        processor_job = wait_for_job(processor_job, ProcessorJob)
        self.assertTrue(processor_job.success)

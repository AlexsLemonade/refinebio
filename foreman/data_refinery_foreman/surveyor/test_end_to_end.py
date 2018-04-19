import time
import json
from django.test import TransactionTestCase, tag
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

def mock_get_sample(accession_code: str):
    """Mock out Sample.objects.get to prevent too much work from being done.

    There are a lot of Samples on this experiment, but we only want to
    do one from each zip file. To accomplish this we make every other
    Sample in the experiment appear to already exist so it is not
    created and processed.
    """
    if accession_code == "GSM1109016" or accession_code == "GSM1108516":
        raise Sample.DoesNotExist

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    # This return value isn't actually used so it doesn't matter what
    # it is. The call to this function is just checking for the
    # existence of the Sample.
    return None


# TransactionTestCase makes database calls complete before the test
# ends.  Otherwise the workers wouldn't actually be able to find the
# job in the database cause it'd be stuck in a transaction.
class ScanUpcEndToEndTestCase(TransactionTestCase):
    @tag("slow")
    @patch('data_refinery_common.models.Sample.objects.get')
    def test_no_op(self, mocked_get_query):
        """Survey, download, then process an experiment we know is NO_OP."""
        # @patch.object(data_refinery_foreman.surveyor.array_express.Sample.objects, "get", mock_get_sample)
        mocked_get_query.side_effect = mock_get_sample
        
        # Prevent a call being made to NCBI's API to determine
        # organism name/id.
        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        survey_job = surveyor.survey_ae_experiment("E-GEOD-45547")

        self.assertTrue(survey_job.success)

        downloader_jobs = DownloaderJob.objects.all()
        logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
        for downloader_job in downloader_jobs:
            downloader_job = wait_for_job(downloader_job, DownloaderJob)
            self.assertTrue(downloader_job.success)

        processor_jobs = ProcessorJob.objects.all()
        logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
        for processor_job in processor_jobs:
            processor_job = wait_for_job(processor_job, ProcessorJob)
            self.assertTrue(processor_job.success)


    # @patch('data_refinery_foreman.surveyor.array_express.requests.get')
    # def test_calls_survey(self, mock_get):
    #     """If source_type is supported calls the appropriate survey method."""
    #     mock_get.side_effect = mocked_requests_get

    #     # Prevent a call being made to NCBI's API to determine
    #     # organism name/id.
    #     organism = Organism(name="HOMO SAPIENS", taxonomy_id=9606, is_scientific_name=True)
    #     organism.save()

    #     survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    #     survey_job.save()
    #     key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
    #                                        key="experiment_accession_code",
    #                                        value="E-GEOD-22166")
    #     key_value_pair.save()

    #     surveyor.run_job(survey_job)
    #     logger.info("Started Survey Job %d, waiting for it to complete.", survey_job.id)
    #     survey_job = wait_for_job(survey_job, SurveyJob)
    #     self.assertTrue(survey_job.success)

    #     batch = Batch.objects.all()[0]
    #     batch = Batch.objects.filter(survey_job=survey_job).get()

    #     downloader_job = batch.downloaderjob_set.get()
    #     logger.info("Survey Job finished, waiting for Downloader Job %d to complete.",
    #                 downloader_job.id)
    #     downloader_job = wait_for_job(downloader_job, DownloaderJob)
    #     self.assertTrue(downloader_job.success)

    #     processor_job = batch.processorjob_set.get()
    #     logger.info("Downloader Job finished, waiting for processor Job %d to complete.",
    #                 processor_job.id)
    #     processor_job = wait_for_job(processor_job, ProcessorJob)
    #     self.assertTrue(processor_job.success)

import glob
import json
import requests
import shutil
import time

from datetime import timedelta, datetime
from django.test import TransactionTestCase, tag
from django.utils import timezone
from unittest.mock import patch, Mock
from test.support import EnvironmentVarGuard # Python >=3

from data_refinery_common.models import (
    Organism,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_foreman.surveyor import surveyor, utils
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment
from data_refinery_foreman.foreman.main import retry_lost_downloader_jobs

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
LOOP_TIME = 5  # seconds
MAX_WAIT_TIME = timedelta(minutes=15)

def wait_for_job(job, job_class: type, start_time: datetime, loop_time: int=None):
    """Monitors the `job_class` table for when `job` is done."""
    loop_time = loop_time if loop_time else LOOP_TIME
    job = job_class.objects.filter(id=job.id).get()
    while job.success is None and timezone.now() - start_time < MAX_WAIT_TIME:
        logger.info("Still polling the %s.",
                    job_class.__name__)
        time.sleep(LOOP_TIME)
        job = job_class.objects.filter(id=job.id).get()

    if timezone.now() - start_time > MAX_WAIT_TIME:
        logger.error("%s job timed out!", job_class.__name__)

    return job


# TransactionTestCase makes database calls complete before the test
# ends.  Otherwise the workers wouldn't actually be able to find the
# job in the database because it'd be stuck in a transaction.
class NoOpEndToEndTestCase(TransactionTestCase):
    @tag("slow")
    def test_no_op(self):
        """Survey, download, then process an experiment we know is NO_OP."""
        # Clear out pre-existing work dirs so there's no conflicts:

        self.env = EnvironmentVarGuard()
        self.env.set('RUNING_IN_CLOUD', 'False')
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Make sure there are no already existing jobs we might poll for unsuccessfully.
            DownloaderJobOriginalFileAssociation.objects.all().delete()
            DownloaderJob.objects.all().delete()
            ProcessorJobOriginalFileAssociation.objects.all().delete()
            ProcessorJob.objects.all().delete()

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            accession_code = "E-GEOD-3303"
            survey_job = surveyor.survey_experiment(accession_code, "ARRAY_EXPRESS")

            self.assertTrue(survey_job.success)

            downloader_jobs = DownloaderJob.objects.all()
            self.assertGreater(downloader_jobs.count(), 0)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            for downloader_job in downloader_jobs:
                downloader_job = wait_for_job(downloader_job, DownloaderJob, start_time)
                self.assertTrue(downloader_job.success)

            processor_jobs = ProcessorJob.objects.all()
            self.assertGreater(processor_jobs.count(), 0)

            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            start_time = timezone.now()
            for processor_job in processor_jobs:
                processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                self.assertTrue(processor_job.success)

            # Test that the unsurveyor deletes all objects related to the experiment
            purge_experiment(accession_code)

            self.assertEqual(Experiment.objects.all().count(), 0)
            self.assertEqual(ExperimentAnnotation.objects.all().count(), 0)
            self.assertEqual(ExperimentSampleAssociation.objects.all().count(), 0)
            self.assertEqual(Sample.objects.all().count(), 0)
            self.assertEqual(SampleAnnotation.objects.all().count(), 0)
            self.assertEqual(OriginalFile.objects.all().count(), 0)
            self.assertEqual(OriginalFileSampleAssociation.objects.all().count(), 0)
            self.assertEqual(SampleResultAssociation.objects.all().count(), 0)
            self.assertEqual(ComputationalResult.objects.all().count(), 0)
            self.assertEqual(ComputationalResultAnnotation.objects.all().count(), 0)
            self.assertEqual(SampleComputedFileAssociation.objects.all().count(), 0)
            self.assertEqual(ComputedFile.objects.all().count(), 0)
            self.assertEqual(DownloaderJob.objects.all().count(), 0)
            self.assertEqual(DownloaderJobOriginalFileAssociation.objects.all().count(), 0)
            self.assertEqual(ProcessorJob.objects.all().count(), 0)
            self.assertEqual(ProcessorJobOriginalFileAssociation.objects.all().count(), 0)


class RedownloadingTestCase(TransactionTestCase):
    @tag("slow")
    def test_array_express_redownloading(self):
        """Survey, download, then process an experiment we know is NO_OP."""
        # Clear out pre-existing work dirs so there's no conflicts:

        self.env = EnvironmentVarGuard()
        self.env.set('RUNING_IN_CLOUD', 'False')
        with self.env:
            for work_dir in glob.glob(LOCAL_ROOT_DIR + "/processor_job_*"):
                shutil.rmtree(work_dir)

            # Make sure there are no already existing jobs we might poll for unsuccessfully.
            DownloaderJobOriginalFileAssociation.objects.all().delete()
            DownloaderJob.objects.all().delete()
            ProcessorJobOriginalFileAssociation.objects.all().delete()
            ProcessorJob.objects.all().delete()

            # Prevent a call being made to NCBI's API to determine
            # organism name/id.
            organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
            organism.save()

            accession_code = "E-GEOD-3303"
            survey_job = surveyor.survey_experiment(accession_code, "ARRAY_EXPRESS")

            self.assertTrue(survey_job.success)

            # This experiment has 12 samples that are contained in the
            # same archive, so only one job is needed.
            downloader_jobs = DownloaderJob.objects.all()
            self.assertEqual(downloader_jobs.count(), 1)

            logger.info("Survey Job finished, waiting for Downloader Jobs to complete.")
            start_time = timezone.now()
            # We want to try and delete the file as quickly as
            # possible, so pass a short loop time and let the waiting
            # loop spin really fast so we lose as little time as
            # possible.
            downloader_job = wait_for_job(downloader_jobs[0], DownloaderJob, start_time, .1)
            self.assertTrue(downloader_job.success)

            # Now we'ere going to delete one of the extracted files but not the other.
            for original_file in OriginalFile.objects.all():
                if not original_file.is_archive:
                    original_file.delete_local_file()
                    break

            # The one downloader job should have extracted 12 files
            # and created 12 processor jobs.
            processor_jobs = ProcessorJob.objects.all()
            self.assertEqual(processor_jobs.count(), 12)

            doomed_processor_job = original_file.processor_jobs.all()[0]
            print(doomed_processor_job.nomad_job_id)
            logger.info(
                "Waiting on processor Nomad job %s to fail because it realized it is missing a file.",
                doomed_processor_job.nomad_job_id
            )

            start_time = timezone.now()
            # This job is somehow succeeeding... LOOK INTO THIS AFTER LUNCH!!!
            with self.assertRaises(ProcessorJob.DoesNotExist):
                wait_for_job(doomed_processor_job, ProcessorJob, start_time)


            # The processor job that had a missing file will have
            # recreated its DownloaderJob, which means there should now be two.
            downloader_jobs = DownloaderJob.objects.all().order_by('-id')
            self.assertEqual(downloader_jobs.count(), 2)

            # However DownloaderJobs don't get queued immediately, so
            # we have to run a foreman function to make it happen:
            retry_lost_downloader_jobs()

            # And we can check that the most recently created
            # DownloaderJob was successful as well:
            recreated_job = downloader_jobs[0]
            recreated_job.refresh_from_db()
            logger.info(
                "Waiting on downloader Nomad job %s",
                recreated_job.nomad_job_id
            )
            recreated_job = wait_for_job(recreated_job, DownloaderJob, start_time)
            self.assertTrue(recreated_job.success)

            # Once the Downloader job succeeds, it should create one
            # and only one processor job, which the total goes back up to 12:
            self.assertEqual(ProcessorJob.objects.all().count(), 12)

            # And finally we can make sure that all 12 of the
            # processor jobs were successful, including the one that
            # got recreated.
            logger.info("Downloader Jobs finished, waiting for processor Jobs to complete.")
            successful_processor_jobs = []
            for processor_job in processor_jobs:
                # One of the two calls to wait_for_job will fail
                # because the job is going to delete itself when it
                # finds that the file it wants to process is missing.
                try:
                    processor_job = wait_for_job(processor_job, ProcessorJob, start_time)
                    if processor_job.success:
                        successful_processor_jobs.append(processor_job)
                except:
                    pass

            # However one of the processor jobs should have succeeded
            # because we didn't delete its file.
            self.assertEqual(len(successful_processor_jobs), 12)

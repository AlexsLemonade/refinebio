import datetime
import math
from test.support import EnvironmentVarGuard  # Python >=3
from unittest.mock import MagicMock, patch

from django.test import TestCase, TransactionTestCase
from django.utils import timezone

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleComputedFileAssociation,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.foreman import job_control

# For use in tests that test the JOB_CREATED_AT_CUTOFF functionality.
DAY_BEFORE_JOB_CUTOFF = job_control.JOB_CREATED_AT_CUTOFF - datetime.timedelta(days=1)

EMPTY_JOB_QUEUE_RESPONSE = {"jobSummaryList": []}


class ForemanTestCase(TestCase):
    def create_downloader_job(self, suffix="e8eaf540"):
        job = DownloaderJob(
            downloader_task="SRA",
            batch_job_id="DOWNLOADER/dispatch-1528945054-" + suffix,
            num_retries=0,
            accession_code="NUNYA",
            success=None,
        )
        job.save()

        og_file = OriginalFile()
        og_file.source_filename = "doesn't matter"
        og_file.filename = "this either"
        og_file.absolute_file_path = "nor this"
        og_file.save()

        assoc1 = DownloaderJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.downloader_job = job
        assoc1.save()

        og_file = OriginalFile()
        og_file.source_filename = "doesn't matter"
        og_file.filename = "this either"
        og_file.absolute_file_path = "nor this"
        og_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = og_file
        assoc.downloader_job = job
        assoc.save()

        return job

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_downloader_job(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_downloader_job()

        job_control.requeue_downloader_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

        self.assertEqual(retried_job.original_files.count(), 2)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_repeated_download_failures(self, mock_list_jobs, mock_send_job):
        """Jobs will be repeatedly retried."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_downloader_job()

        for i in range(job_control.MAX_NUM_RETRIES):
            job_control.handle_downloader_jobs([job])
            self.assertEqual(i + 1, len(mock_send_job.mock_calls))

            jobs = DownloaderJob.objects.all().order_by("-id")
            previous_job = jobs[1]
            self.assertTrue(previous_job.retried)
            self.assertEqual(previous_job.num_retries, i)
            self.assertFalse(previous_job.success)

            job = jobs[0]
            self.assertFalse(job.retried)
            self.assertEqual(job.num_retries, i + 1)

        # Once MAX_NUM_RETRIES has been hit handle_repeated_failure
        # should be called.
        job_control.handle_downloader_jobs([job])
        last_job = DownloaderJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, job_control.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_retrying_failed_downloader_jobs(self, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_downloader_job()
        job.success = False
        job.save()

        job_control.retry_failed_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_hung_downloader_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "FAILED"}]
        }

        job = self.create_downloader_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        # This second job is no longer listed by Batch, needs to be retried.
        job2 = self.create_downloader_job()
        job2.start_time = timezone.now()
        job2.batch_job_id = "MISSING"
        job2.save()

        job_control.retry_hung_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_hung_downloader_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Tests that we don't restart downloader jobs that are still running."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNING"}]
        }

        job = self.create_downloader_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_hung_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_lost_downloader_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_downloader_job()
        job.created_at = timezone.now()
        job.save()

        # Test that jobs with batch_job_id still get requeued.
        job2 = self.create_downloader_job()
        job2.created_at = timezone.now()
        job2.batch_job_id = "MISSING"
        job2.save()

        job_control.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_old_downloader_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_downloader_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        job_control.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        self.assertEqual(1, DownloaderJob.objects.all().count())

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_lost_downloader_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Make sure that we don't retry downloader jobs we shouldn't."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNABLE"}]
        }

        job = self.create_downloader_job()
        job.created_at = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    def create_processor_job(self, pipeline="AFFY_TO_PCL", ram_amount=2048, start_time=None):
        job = ProcessorJob(
            pipeline_applied=pipeline,
            batch_job_id="PROCESSOR/dispatch-1528945054-e8eaf540",
            ram_amount=ram_amount,
            num_retries=0,
            volume_index="1",
            success=None,
            start_time=start_time,
        )
        job.save()

        og_file = OriginalFile()
        og_file.source_filename = "doesn't matter"
        og_file.filename = "this either"
        og_file.absolute_file_path = "nor this"
        og_file.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        og_file = OriginalFile()
        og_file.source_filename = "doesn't matter"
        og_file.filename = "this either"
        og_file.absolute_file_path = "nor this"
        og_file.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = og_file
        assoc.processor_job = job
        assoc.save()

        return job

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_processor_job(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job()

        job_control.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_processor_job_no_volume(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job()
        job.volume_index = None
        job.save()

        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "True")

        with self.settings(RUNNING_IN_CLOUD=True):
            job_control.requeue_processor_job(job)

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertEqual(retried_job.volume_index, "0")

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_compendia_job_no_volume(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job()
        job.volume_index = None
        job.pipeline_applied = "CREATE_COMPENDIA"
        job.save()

        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "True")

        with self.settings(RUNNING_IN_CLOUD=True):
            job_control.requeue_processor_job(job)

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertEqual(retried_job.volume_index, None)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_processor_job_w_more_ram(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job(
            pipeline="SALMON", ram_amount=16384, start_time=timezone.now()
        )

        job_control.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)
        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertEqual(original_job.ram_amount, 16384)
        self.assertEqual(retried_job.ram_amount, 32768)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_repeated_processor_failures(self, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        """Jobs will be repeatedly retried."""
        job = self.create_processor_job()

        for i in range(job_control.MAX_NUM_RETRIES):
            job_control.handle_processor_jobs([job])
            self.assertEqual(i + 1, len(mock_send_job.mock_calls))

            jobs = ProcessorJob.objects.all().order_by("-id")
            previous_job = jobs[1]
            self.assertTrue(previous_job.retried)
            self.assertEqual(previous_job.num_retries, i)
            self.assertFalse(previous_job.success)

            job = jobs[0]
            self.assertFalse(job.retried)
            self.assertEqual(job.num_retries, i + 1)

        # Once MAX_NUM_RETRIES has been hit handle_repeated_failure
        # should be called.
        job_control.handle_processor_jobs([job])
        last_job = ProcessorJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, job_control.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_retrying_failed_processor_jobs(self, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_processor_job()
        job.success = False
        job.save()

        job_control.retry_failed_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_hung_processor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "FAILED"}]
        }

        job = self.create_processor_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job2 = self.create_processor_job()
        job2.start_time = timezone.now()
        job2.batch_job_id = "MISSING"
        job2.save()

        job_control.retry_hung_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_hung_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Tests that we don't restart processor jobs that are still running."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNING"}]
        }

        job = self.create_processor_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_hung_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_lost_processor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_processor_job()
        job.save()

        job2 = self.create_processor_job()
        job2.batch_job_id = "MISSING"
        job2.save()

        job_control.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_lost_smasher_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        """Make sure that the smasher jobs will get retried even though they
        don't have a volume_index.

        I'm not entirely sure this test is still necessary but we'll
        need a separate smasher compute environment so this could test
        that once it's done.
        """
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_processor_job(pipeline="SMASHER")
        job.volume_index = None  # Smasher jobs won't have a volume_index.
        job.save()

        job_control.retry_lost_processor_jobs()

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_old_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_processor_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        job_control.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        self.assertEqual(1, ProcessorJob.objects.all().count())

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_lost_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Make sure that we don't retry processor jobs we shouldn't."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNABLE"}]
        }

        job = self.create_processor_job()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_janitor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_processor_job(pipeline="JANITOR")
        job.save()

        job_control.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        self.assertEqual(len(jobs), 1)

    def create_survey_job(self):
        job = SurveyJob(
            source_type="SRA",
            batch_job_id="SURVEYOR/dispatch-1528945054-e8eaf540",
            num_retries=0,
            success=None,
        )

        job.save()

        sjkv = SurveyJobKeyValue()
        sjkv.key = "experiment_accession_code"
        sjkv.value = "RJ-1234-XYZ"
        sjkv.survey_job = job
        sjkv.save()

        return job

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_requeuing_survey_job(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_survey_job()

        job_control.requeue_survey_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_repeated_survey_failures(self, mock_list_jobs, mock_send_job):
        """Jobs will be repeatedly retried."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_survey_job()

        for i in range(job_control.MAX_NUM_RETRIES):
            job_control.handle_survey_jobs([job])
            self.assertEqual(i + 1, len(mock_send_job.mock_calls))

            jobs = SurveyJob.objects.all().order_by("-id")
            previous_job = jobs[1]
            self.assertTrue(previous_job.retried)
            self.assertEqual(previous_job.num_retries, i)
            self.assertFalse(previous_job.success)

            job = jobs[0]
            self.assertFalse(job.retried)
            self.assertEqual(job.num_retries, i + 1)

        # Once MAX_NUM_RETRIES has been hit handle_repeated_failure
        # should be called.
        job_control.handle_survey_jobs([job])
        last_job = SurveyJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, job_control.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

        # MAX TOTAL tests
        self.env = EnvironmentVarGuard()
        self.env.set("MAX_TOTAL_JOBS", "0")
        with self.env:
            job = self.create_survey_job()
            result = job_control.handle_survey_jobs([job])
            self.assertFalse(result)

        self.env.set("MAX_TOTAL_JOBS", "1000")
        with self.env:
            job = self.create_survey_job()
            result = job_control.requeue_survey_job(job)
            self.assertTrue(result)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    def test_retrying_failed_survey_jobs(self, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_survey_job()
        job.success = False
        job.save()

        job_control.retry_failed_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_hung_survey_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_survey_job()
        job.start_time = timezone.now()
        job.save()

        job2 = self.create_survey_job()
        job2.start_time = timezone.now()
        job2.batch_job_id = "MISSING"
        job2.save()

        job_control.retry_hung_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        jobs = SurveyJob.objects.order_by("id")
        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_hung_survey_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        """Tests that we don't restart survey jobs that are still running."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNING"}]
        }

        job = self.create_survey_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_hung_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_retrying_lost_survey_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_survey_job()
        job.created_at = timezone.now()
        job.save()

        job2 = self.create_survey_job()
        job2.created_at = timezone.now()
        job2.save()

        job_control.retry_lost_survey_jobs()

        self.assertEqual(len(mock_send_job.mock_calls), 2)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        jobs = SurveyJob.objects.order_by("id")
        original_job2 = jobs[1]
        self.assertTrue(original_job2.retried)
        self.assertEqual(original_job2.num_retries, 0)
        self.assertFalse(original_job2.success)

        retried_job = jobs[2]
        self.assertEqual(retried_job.num_retries, 1)

        retried_job2 = jobs[3]
        self.assertEqual(retried_job2.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_old_survey_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE

        job = self.create_survey_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        job_control.retry_lost_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        self.assertEqual(1, SurveyJob.objects.all().count())

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    @patch("data_refinery_foreman.foreman.job_control.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.job_control.batch.describe_jobs")
    def test_not_retrying_lost_survey_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        """Make sure that we don't retry survey jobs we shouldn't."""
        mock_send_job.return_value = True
        mock_list_jobs.return_value = EMPTY_JOB_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {
            "jobSummaryList": [{"jobId": "FINDME", "status": "RUNNABLE"}]
        }

        job = self.create_survey_job()
        job.created_at = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job_control.retry_lost_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_control.send_job")
    def test_janitor(self, mock_send_job):
        """For now super simple, this will evolve once we have more than one compute environment"""
        mock_send_job.return_value = True

        job_control.send_janitor_jobs()

        self.assertEqual(ProcessorJob.objects.all().count(), 1)


class CleanDatabaseTestCase(TransactionTestCase):
    def test_cleandb(self):
        sample = Sample()
        sample.save()

        result = ComputationalResult()
        result.save()

        good_file = ComputedFile()
        good_file.s3_bucket = "my_cool_bucket"
        good_file.s3_key = "my_sweet_key"
        good_file.size_in_bytes = 1337
        good_file.result = result
        good_file.is_public = True
        good_file.is_smashable = True
        good_file.save()

        sca = SampleComputedFileAssociation()
        sca.sample = sample
        sca.computed_file = good_file
        sca.save()

        bad_file = ComputedFile()
        bad_file.s3_bucket = None
        bad_file.s3_key = None
        bad_file.result = result
        bad_file.size_in_bytes = 7331
        bad_file.is_public = True
        bad_file.is_smashable = True
        bad_file.save()

        sca = SampleComputedFileAssociation()
        sca.sample = sample
        sca.computed_file = bad_file
        sca.save()

        self.assertEqual(sample.computed_files.count(), 2)
        self.assertEqual(sample.get_most_recent_smashable_result_file().id, bad_file.id)
        job_control.clean_database()
        self.assertEqual(sample.get_most_recent_smashable_result_file().id, good_file.id)

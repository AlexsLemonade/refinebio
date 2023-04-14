from test.support import EnvironmentVarGuard
from unittest.mock import patch

from django.conf import settings
from django.test import TestCase
from django.utils import timezone

from tests.foreman.utils import create_downloader_job, create_processor_job, create_survey_job

from data_refinery_common.models import DownloaderJob, ProcessorJob, SurveyJob
from data_refinery_foreman.foreman import job_requeuing


def fake_send_job(job_type, job, is_dispatch=False):
    job.batch_job_queue = settings.AWS_BATCH_QUEUE_WORKERS_NAMES[0]
    job.save()

    return True


class JobRequeuingTestCase(TestCase):
    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_downloader_job(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_downloader_job()

        job_requeuing.requeue_downloader_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

        self.assertEqual(retried_job.original_files.count(), 2)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_processor_job(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_processor_job()

        job_requeuing.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_processor_job_no_batch_job_queue(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_processor_job()
        job.batch_job_queue = None
        job.save()

        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "True")

        with self.settings(RUNNING_IN_CLOUD=True):
            job_requeuing.requeue_processor_job(job)

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertIsNotNone(retried_job.batch_job_queue)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_compendia_job_no_batch_job_queue(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_processor_job()
        job.batch_job_queue = None
        job.pipeline_applied = "CREATE_COMPENDIA"
        job.save()

        self.env = EnvironmentVarGuard()
        self.env.set("RUNING_IN_CLOUD", "True")

        with self.settings(RUNNING_IN_CLOUD=True):
            job_requeuing.requeue_processor_job(job)

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertIsNotNone(retried_job.batch_job_queue)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_processor_job_w_more_ram(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_processor_job(pipeline="SALMON", ram_amount=16384, start_time=timezone.now())

        job_requeuing.requeue_processor_job(job)
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

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    def test_requeuing_survey_job(self, mock_send_job):
        mock_send_job.side_effect = fake_send_job

        job = create_survey_job()

        job_requeuing.requeue_survey_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

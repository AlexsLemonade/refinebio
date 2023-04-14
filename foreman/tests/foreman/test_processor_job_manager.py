import datetime
from unittest.mock import patch

from django.conf import settings
from django.test import TestCase
from django.utils import timezone

from tests.foreman.utils import create_processor_job

from data_refinery_common.models import ProcessorJob
from data_refinery_foreman.foreman import processor_job_manager, utils

# For use in tests that test the JOB_CREATED_AT_CUTOFF functionality.
DAY_BEFORE_JOB_CUTOFF = utils.JOB_CREATED_AT_CUTOFF - datetime.timedelta(days=1)

EMPTY_LIST_JOBS_QUEUE_RESPONSE = {"jobSummaryList": []}
EMPTY_DESCRIBE_JOBS_QUEUE_RESPONSE = {"jobs": []}


def fake_send_job(job_type, job, is_dispatch=False):
    job.batch_job_queue = settings.AWS_BATCH_QUEUE_WORKERS_NAMES[0]
    job.save()

    return True


class ProcessorJobManagerTestCase(TestCase):
    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    def test_repeated_processor_failures(self, mock_list_jobs, mock_send_job):
        """Jobs will be repeatedly retried."""
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE

        job = create_processor_job()

        for i in range(utils.MAX_NUM_RETRIES):
            processor_job_manager.handle_processor_jobs([job])
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
        processor_job_manager.handle_processor_jobs([job])
        last_job = ProcessorJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, utils.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    def test_retrying_failed_processor_jobs(self, mock_list_jobs, mock_send_job):
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE

        job = create_processor_job()
        job.success = False
        job.save()

        processor_job_manager.retry_failed_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_retrying_hung_processor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {"jobs": [{"jobId": "FINDME", "status": "FAILED"}]}

        job = create_processor_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        job2 = create_processor_job()
        job2.start_time = timezone.now()
        job2.batch_job_id = "MISSING"
        job2.save()

        processor_job_manager.retry_hung_processor_jobs()
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

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_not_retrying_hung_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Tests that we don't restart processor jobs that are still running."""
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {"jobs": [{"jobId": "FINDME", "status": "RUNNING"}]}

        job = create_processor_job()
        job.start_time = timezone.now()
        job.batch_job_id = "FINDME"
        job.save()

        processor_job_manager.retry_hung_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_retrying_lost_processor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_DESCRIBE_JOBS_QUEUE_RESPONSE

        job = create_processor_job()
        job.save()

        job2 = create_processor_job()
        job2.batch_job_id = "MISSING"
        job2.save()

        processor_job_manager.retry_lost_processor_jobs()
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

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_retrying_lost_smasher_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        """Make sure that the smasher jobs will get retried even though they
        don't have a volume_index.

        I'm not entirely sure this test is still necessary but we'll
        need a separate smasher compute environment so this could test
        that once it's done.
        """
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_DESCRIBE_JOBS_QUEUE_RESPONSE

        job = create_processor_job(pipeline="SMASHER")
        job.volume_index = None  # Smasher jobs won't have a volume_index.
        job.save()

        processor_job_manager.retry_lost_processor_jobs()

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_not_retrying_old_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_DESCRIBE_JOBS_QUEUE_RESPONSE

        job = create_processor_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        processor_job_manager.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        self.assertEqual(1, ProcessorJob.objects.all().count())

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_not_retrying_lost_processor_jobs(
        self, mock_describe_jobs, mock_list_jobs, mock_send_job
    ):
        """Make sure that we don't retry processor jobs we shouldn't."""
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = {"jobs": [{"jobId": "FINDME", "status": "RUNNABLE"}]}

        job = create_processor_job()
        job.batch_job_id = "FINDME"
        job.save()

        processor_job_manager.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch("data_refinery_foreman.foreman.job_requeuing.send_job")
    @patch("data_refinery_common.message_queue.batch.list_jobs")
    @patch("data_refinery_foreman.foreman.utils.batch.describe_jobs")
    def test_not_retrying_janitor_jobs(self, mock_describe_jobs, mock_list_jobs, mock_send_job):
        mock_send_job.side_effect = fake_send_job
        mock_list_jobs.return_value = EMPTY_LIST_JOBS_QUEUE_RESPONSE
        mock_describe_jobs.return_value = EMPTY_DESCRIBE_JOBS_QUEUE_RESPONSE

        job = create_processor_job(pipeline="JANITOR")
        job.save()

        processor_job_manager.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by("id")
        self.assertEqual(len(jobs), 1)

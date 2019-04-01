from unittest.mock import patch, MagicMock
import datetime
import time
from django.utils import timezone
from django.test import TestCase
from data_refinery_foreman.foreman import main
from data_refinery_common.models import (
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
    SurveyJob,
    SurveyJobKeyValue,
)
from test.support import EnvironmentVarGuard # Python >=3

# For use in tests that test the JOB_CREATED_AT_CUTOFF functionality.
DAY_BEFORE_JOB_CUTOFF = main.JOB_CREATED_AT_CUTOFF - datetime.timedelta(days=1)

class ForemanTestCase(TestCase):
    def create_downloader_job(self, suffix="e8eaf540"):
        job = DownloaderJob(downloader_task="SRA",
                            nomad_job_id="DOWNLOADER/dispatch-1528945054-" + suffix,
                            num_retries=0,
                            accession_code="NUNYA",
                            success=None)
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

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_downloader_job(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_downloader_job()

        main.requeue_downloader_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

        self.assertEqual(retried_job.original_files.count(), 2)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_repeated_download_failures(self, mock_send_job):
        """Jobs will be repeatedly retried."""
        mock_send_job.return_value = True

        job = self.create_downloader_job()

        for i in range(main.MAX_NUM_RETRIES):
            main.handle_downloader_jobs([job])
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
        main.handle_downloader_jobs([job])
        last_job = DownloaderJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, main.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_failed_downloader_jobs(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_downloader_job()
        job.success = False
        job.save()

        main.retry_failed_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_many_failed_downloader_jobs(self, mock_send_job):
        mock_send_job.return_value = True

        for x in range(0, 10000):
            job = self.create_downloader_job(str(x))
            job.success = False
            job.save()

        main.retry_failed_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 10000)

        jobs = DownloaderJob.objects.order_by('id')

        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[10001]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_hung_downloader_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_downloader_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_hung_downloader_jobs(self, mock_nomad, mock_send_job):
        """Tests that we don't restart downloader jobs that are still running."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "running"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_downloader_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_lost_downloader_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_downloader_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_old_downloader_jobs(self, mock_nomad, mock_send_job):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_downloader_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        main.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = DownloaderJob.objects.order_by('id')
        self.assertEqual(1, DownloaderJob.objects.all().count())

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_lost_downloader_jobs_time(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_downloader_job()
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + datetime.timedelta(minutes=1))
        job.save()

        main.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_lost_downloader_jobs(self, mock_nomad, mock_send_job):
        """Make sure that we don't retry downloader jobs we shouldn't."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "pending"}
            return ret_value

        mock_nomad.side_effect=mock_init_nomad

        job = self.create_downloader_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_downloader_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    def create_processor_job(self, pipeline="AFFY_TO_PCL", ram_amount=2048):
        job = ProcessorJob(pipeline_applied=pipeline,
                           nomad_job_id="PROCESSOR/dispatch-1528945054-e8eaf540",
                           ram_amount=ram_amount,
                           num_retries=0,
                           volume_index="1",
                           success=None)
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

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_processor_job(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        job = self.create_processor_job()

        main.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_processor_job_w_more_ram(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        job = self.create_processor_job(pipeline="SALMON", ram_amount=16384)

        main.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)
        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertEqual(original_job.ram_amount, 16384)
        self.assertEqual(retried_job.ram_amount, 32768)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_repeated_processor_failures(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        """Jobs will be repeatedly retried."""
        job = self.create_processor_job()

        for i in range(main.MAX_NUM_RETRIES):
            main.handle_processor_jobs([job])
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
        main.handle_processor_jobs([job])
        last_job = ProcessorJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, main.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_failed_processor_jobs(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        job = self.create_processor_job()
        job.success = False
        job.save()

        main.retry_failed_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_not_retrying_wrong_volume_index(self, mock_send_job, mock_get_active_volumes):
        """If a volume isn't mounted then we shouldn't queue jobs for it."""
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"2", "3"}

        job = self.create_processor_job()
        job.success = False
        job.save()

        main.retry_failed_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        self.assertEqual(len(jobs), 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_hung_processor_jobs(self, mock_nomad, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_processor_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_hung_processor_jobs(self, mock_nomad, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        """Tests that we don't restart processor jobs that are still running."""
        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "running"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_processor_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_lost_processor_jobs(self, mock_nomad, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_processor_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_processor_jobs()

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_old_processor_jobs(self, mock_nomad, mock_send_job):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_processor_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        main.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        self.assertEqual(1, ProcessorJob.objects.all().count())

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_lost_processor_jobs(self, mock_nomad, mock_send_job, mock_get_active_volumes):
        """Make sure that we don't retry processor jobs we shouldn't."""
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "pending"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_processor_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_lost_processor_jobs_time(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        job = self.create_processor_job()
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + datetime.timedelta(minutes=1))
        job.save()

        main.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)


    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_not_retrying_janitor_jobs(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        job = self.create_processor_job(pipeline="JANITOR")
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + datetime.timedelta(minutes=1))
        job.save()

        main.retry_lost_processor_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        self.assertEqual(len(jobs), 1)

    def create_survey_job(self):
        job = SurveyJob(source_type="SRA",
                        nomad_job_id="SURVEYOR/dispatch-1528945054-e8eaf540",
                        num_retries=0, success=None)

        job.save()

        sjkv = SurveyJobKeyValue()
        sjkv.key = "experiment_accession_code"
        sjkv.value = "RJ-1234-XYZ"
        sjkv.survey_job = job
        sjkv.save()

        return job

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_survey_job(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_survey_job()

        main.requeue_survey_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_repeated_survey_failures(self, mock_send_job):
        """Jobs will be repeatedly retried."""
        mock_send_job.return_value = True

        job = self.create_survey_job()

        for i in range(main.MAX_NUM_RETRIES):
            main.handle_survey_jobs([job])
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
        main.handle_survey_jobs([job])
        last_job = SurveyJob.objects.all().order_by("-id")[0]
        self.assertTrue(last_job.retried)
        self.assertEqual(last_job.num_retries, main.MAX_NUM_RETRIES)
        self.assertFalse(last_job.success)

        # MAX TOTAL tests
        self.env = EnvironmentVarGuard()
        self.env.set('MAX_TOTAL_JOBS', '0')
        with self.env:
            job = self.create_survey_job()
            result = main.handle_survey_jobs([job])
            self.assertFalse(result)

        self.env.set('MAX_TOTAL_JOBS', '1000')
        with self.env:
            job = self.create_survey_job()
            result = main.requeue_survey_job(job)
            self.assertTrue(result)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_failed_survey_jobs(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_survey_job()
        job.success = False
        job.save()

        main.retry_failed_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_hung_survey_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_survey_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_hung_survey_jobs(self, mock_nomad, mock_send_job):
        """Tests that we don't restart survey jobs that are still running."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "running"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_survey_job()
        job.start_time = timezone.now()
        job.save()

        main.retry_hung_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_lost_survey_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_survey_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_survey_jobs()

        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_old_survey_jobs(self, mock_nomad, mock_send_job):
        """Makes sure temporary logic to limit the Foreman's scope works."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "dead"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_survey_job()
        job.created_at = DAY_BEFORE_JOB_CUTOFF
        job.save()

        main.retry_lost_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = SurveyJob.objects.order_by('id')
        self.assertEqual(1, SurveyJob.objects.all().count())

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_not_retrying_lost_survey_jobs(self, mock_nomad, mock_send_job):
        """Make sure that we don't retry survey jobs we shouldn't."""
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = lambda _: {"Status": "pending"}
            return ret_value

        mock_nomad.side_effect = mock_init_nomad

        job = self.create_survey_job()
        job.created_at = timezone.now()
        job.save()

        main.retry_lost_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_lost_survey_jobs_time(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_survey_job()
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + datetime.timedelta(minutes=1))
        job.save()

        main.retry_lost_survey_jobs()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.get_active_volumes')
    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_janitor(self, mock_send_job, mock_get_active_volumes):
        mock_send_job.return_value = True
        mock_get_active_volumes.return_value = {"1", "2", "3"}

        for p in ["1", "2", "3"]:
            pj = ProcessorJob()
            pj.volume_index = p
            pj.save()

        main.send_janitor_jobs()

        self.assertEqual(ProcessorJob.objects.all().count(), 6)
        self.assertEqual(ProcessorJob.objects.filter(pipeline_applied="JANITOR").count(), 3)

        # Make sure that the janitors are dispatched to the correct volumes.
        ixs = ["1", "2", "3"]
        for p in ProcessorJob.objects.filter(pipeline_applied="JANITOR"):
            self.assertTrue(p.volume_index in ixs)
            ixs.remove(p.volume_index)

    def test_get_max_downloader_jobs(self):
        self.assertNotEqual(main.get_max_downloader_jobs(), 0)

# class JobPrioritizationTestCase(TestCase):
#     def setUp(self):
#         """Create a lot of resources that could be associated with either
#         ProcessorJobs or DownloaderJobs. Since the logic of when to actually
#         queue these is the same, we can use these for testing both. However
#         The actual jobs that will be queued need to be created by the job-type
#         specific functions.
#         """
#         human = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
#         human.save()
#         zebrafish = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
#         zebrafish.save()

#         # Salmon experiment that is 50% complete.
#         experiment = Experiment(accession_code='ERP036000')
#         experiment.save()

#         ## First sample, this one has been processed.
#         pj = ProcessorJob()
#         pj.accession_code = "ERR036000"
#         pj.pipeline_applied = "SALMON"
#         pj.success = True
#         pj.save()

#         og = OriginalFile()
#         og.filename = "ERR036000.fastq.gz"
#         og.source_filename = "ERR036000.fastq.gz"
#         og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz"
#         og.is_archive = True
#         og.save()

#         sample = Sample()
#         sample.accession_code = 'ERR036000'
#         sample.organism = human
#         sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = sample
#         assoc.original_file = og
#         assoc.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = pj
#         assoc.original_file = og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = sample
#         assoc.experiment = experiment
#         assoc.save()

#         ## Second sample, this one hasn't been processed.
#         self.in_progress_salmon_og = OriginalFile()
#         self.in_progress_salmon_og.filename = "ERR036001.fastq.gz"
#         self.in_progress_salmon_og.source_filename = "ERR036001.fastq.gz"
#         self.in_progress_salmon_og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036001/ERR036001_1.fastq.gz"
#         self.in_progress_salmon_og.is_archive = True
#         self.in_progress_salmon_og.save()

#         self.in_progress_salmon_sample = Sample()
#         self.in_progress_salmon_sample.accession_code = 'ERR036001'
#         self.in_progress_salmon_sample.organism = human
#         self.in_progress_salmon_sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = self.in_progress_salmon_sample
#         assoc.original_file = self.in_progress_salmon_og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = self.in_progress_salmon_sample
#         assoc.experiment = experiment
#         assoc.save()


#         # Salmon experiment that is 0% complete.
#         experiment = Experiment(accession_code='ERP037000')
#         experiment.save()

#         self.unstarted_salmon_og = OriginalFile()
#         self.unstarted_salmon_og.filename = "ERR037001.fastq.gz"
#         self.unstarted_salmon_og.source_filename = "ERR037001.fastq.gz"
#         self.unstarted_salmon_og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR037/ERR037001/ERR037001_1.fastq.gz"
#         self.unstarted_salmon_og.is_archive = True
#         self.unstarted_salmon_og.save()

#         self.unstarted_salmon_sample = Sample()
#         self.unstarted_salmon_sample.accession_code = 'ERR037001'
#         self.unstarted_salmon_sample.organism = human
#         self.unstarted_salmon_sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = self.unstarted_salmon_sample
#         assoc.original_file = self.unstarted_salmon_og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = self.unstarted_salmon_sample
#         assoc.experiment = experiment
#         assoc.save()


#         # Zebrafish experiment.
#         experiment = Experiment(accession_code='ERP038000')
#         experiment.save()

#         self.zebrafish_og = OriginalFile()
#         self.zebrafish_og.source_filename = "ERR038001.fastq.gz"
#         self.zebrafish_og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR038/ERR038001/ERR038001_1.fastq.gz"
#         self.zebrafish_og.is_archive = True
#         self.zebrafish_og.save()

#         self.zebrafish_sample = Sample()
#         self.zebrafish_sample.accession_code = 'ERR038001'
#         self.zebrafish_sample.organism = zebrafish
#         self.zebrafish_sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = self.zebrafish_sample
#         assoc.original_file = self.zebrafish_og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = self.zebrafish_sample
#         assoc.experiment = experiment
#         assoc.save()


#         # Pediatric experiment.
#         experiment = Experiment(accession_code='GSE100568')
#         experiment.save()

#         self.pediatric_og = OriginalFile()
#         self.pediatric_og.source_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100568&format=file"
#         self.pediatric_og.is_archive = True
#         self.pediatric_og.save()

#         self.pediatric_sample = Sample()
#         self.pediatric_sample.accession_code = 'GSM2687180'
#         self.pediatric_sample.organism = human
#         self.pediatric_sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = self.pediatric_sample
#         assoc.original_file = self.pediatric_og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = self.pediatric_sample
#         assoc.experiment = experiment
#         assoc.save()


#         # hgu133plus2 experiment.
#         experiment = Experiment(accession_code='GSE100014')
#         experiment.save()

#         self.hgu133plus2_og = OriginalFile()
#         self.hgu133plus2_og.source_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100014&format=file"
#         self.hgu133plus2_og.is_archive = True
#         self.hgu133plus2_og.save()

#         self.hgu133plus2_sample = Sample()
#         self.hgu133plus2_sample.accession_code = 'GSM2667926'
#         self.hgu133plus2_sample.organism = human
#         self.hgu133plus2_sample.save()

#         assoc = OriginalFileSampleAssociation()
#         assoc.sample = self.hgu133plus2_sample
#         assoc.original_file = self.hgu133plus2_og
#         assoc.save()

#         assoc = ExperimentSampleAssociation()
#         assoc.sample = self.hgu133plus2_sample
#         assoc.experiment = experiment
#         assoc.save()

#     @patch('data_refinery_foreman.foreman.main.Nomad')
#     @patch('data_refinery_foreman.foreman.main.requeue_downloader_job')
#     def test_handle_downloader_jobs(self, mock_requeue_downloader_job, mock_nomad):
#         """Tests the prioritization of downloader jobs.

#         We want zebrafish jobs to be first, then jobs for hgu133plus2,
#         then jobs for pediatric cancer, finally salmon jobs should be
#         prioritized based on how close to completion they are."""

#         def mock_init_nomad(host, port=0, timeout=0):
#             ret_value = MagicMock()
#             ret_value.jobs = MagicMock()
#             ret_value.jobs.get_jobs = MagicMock()
#             ret_value.jobs.get_jobs.side_effect = lambda: []
#             return ret_value

#         mock_nomad.side_effect = mock_init_nomad

#         unstarted_salmon_job = DownloaderJob()
#         unstarted_salmon_job.accession_code = self.unstarted_salmon_sample.accession_code
#         unstarted_salmon_job.save()

#         assoc = DownloaderJobOriginalFileAssociation()
#         assoc.downloader_job = unstarted_salmon_job
#         assoc.original_file = self.unstarted_salmon_og
#         assoc.save()

#         in_progress_salmon_job = DownloaderJob()
#         in_progress_salmon_job.accession_code = self.in_progress_salmon_sample.accession_code
#         in_progress_salmon_job.save()

#         assoc = DownloaderJobOriginalFileAssociation()
#         assoc.downloader_job = in_progress_salmon_job
#         assoc.original_file = self.in_progress_salmon_og
#         assoc.save()

#         zebrafish_job = DownloaderJob()
#         zebrafish_job.accession_code = self.zebrafish_sample.accession_code
#         zebrafish_job.save()

#         assoc = DownloaderJobOriginalFileAssociation()
#         assoc.downloader_job = zebrafish_job
#         assoc.original_file = self.zebrafish_og
#         assoc.save()

#         pediatric_job = DownloaderJob()
#         pediatric_job.accession_code = self.pediatric_sample.accession_code
#         pediatric_job.save()

#         assoc = DownloaderJobOriginalFileAssociation()
#         assoc.downloader_job = pediatric_job
#         assoc.original_file = self.pediatric_og
#         assoc.save()

#         hgu133plus2_job = DownloaderJob()
#         hgu133plus2_job.accession_code = self.hgu133plus2_sample.accession_code
#         hgu133plus2_job.save()

#         assoc = DownloaderJobOriginalFileAssociation()
#         assoc.downloader_job = hgu133plus2_job
#         assoc.original_file = self.hgu133plus2_og
#         assoc.save()

#         jobs = [unstarted_salmon_job,
#                 in_progress_salmon_job,
#                 hgu133plus2_job,
#                 zebrafish_job,
#                 pediatric_job
#         ]
#         jobs_in_correct_order = [zebrafish_job,
#                                  hgu133plus2_job,
#                                  pediatric_job,
#                                  in_progress_salmon_job,
#                                  unstarted_salmon_job
#         ]

#         main.handle_downloader_jobs(jobs)

#         for count, job in enumerate(jobs_in_correct_order):
#             # Calls are a weird object that I think is just basically
#             # a tuple. Index 1 of a call object is the arguments
#             # tuple, we're interested in the first argument
#             job_called_at_count = mock_requeue_downloader_job.mock_calls[count][1][0]
#             self.assertEqual(job.id, job_called_at_count.id)


#     @patch('data_refinery_foreman.foreman.main.Nomad')
#     @patch('data_refinery_foreman.foreman.main.requeue_processor_job')
#     def test_handle_processor_jobs(self, mock_requeue_processor_job, mock_nomad):
#         """Tests the prioritization of processor jobs.

#         We want zebrafish jobs to be first, then jobs for hgu133plus2,
#         then jobs for pediatric cancer, finally salmon jobs should be
#         prioritized based on how close to completion they are."""

#         def mock_init_nomad(host, port=0, timeout=0):
#             ret_value = MagicMock()
#             ret_value.jobs = MagicMock()
#             ret_value.jobs.get_jobs = MagicMock()
#             ret_value.jobs.get_jobs.side_effect = lambda: []
#             return ret_value

#         mock_nomad.side_effect = mock_init_nomad

#         unstarted_salmon_job = ProcessorJob()
#         unstarted_salmon_job.accession_code = self.unstarted_salmon_sample.accession_code
#         unstarted_salmon_job.pipeline_applied = "SALMON"
#         unstarted_salmon_job.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = unstarted_salmon_job
#         assoc.original_file = self.unstarted_salmon_og
#         assoc.save()

#         in_progress_salmon_job = ProcessorJob()
#         in_progress_salmon_job.accession_code = self.in_progress_salmon_sample.accession_code
#         in_progress_salmon_job.pipeline_applied = "SALMON"
#         in_progress_salmon_job.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = in_progress_salmon_job
#         assoc.original_file = self.in_progress_salmon_og
#         assoc.save()

#         zebrafish_job = ProcessorJob()
#         zebrafish_job.accession_code = self.zebrafish_sample.accession_code
#         zebrafish_job.pipeline_applied = "SALMON"
#         zebrafish_job.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = zebrafish_job
#         assoc.original_file = self.zebrafish_og
#         assoc.save()

#         pediatric_job = ProcessorJob()
#         pediatric_job.accession_code = self.pediatric_sample.accession_code
#         pediatric_job.pipeline_applied = "SALMON"
#         pediatric_job.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = pediatric_job
#         assoc.original_file = self.pediatric_og
#         assoc.save()

#         hgu133plus2_job = ProcessorJob()
#         hgu133plus2_job.accession_code = self.hgu133plus2_sample.accession_code
#         hgu133plus2_job.pipeline_applied = "SALMON"
#         hgu133plus2_job.save()

#         assoc = ProcessorJobOriginalFileAssociation()
#         assoc.processor_job = hgu133plus2_job
#         assoc.original_file = self.hgu133plus2_og
#         assoc.save()

#         jobs = [unstarted_salmon_job,
#                 in_progress_salmon_job,
#                 hgu133plus2_job,
#                 zebrafish_job,
#                 pediatric_job
#         ]
#         jobs_in_correct_order = [zebrafish_job,
#                                  hgu133plus2_job,
#                                  pediatric_job,
#                                  in_progress_salmon_job,
#                                  unstarted_salmon_job
#         ]

#         main.handle_processor_jobs(jobs)

#         for count, job in enumerate(jobs_in_correct_order):
#             # Calls are a weird object that I think is just basically
#             # a tuple. Index 1 of a call object is the arguments
#             # tuple, we're interested in the first argument
#             job_called_at_count = mock_requeue_processor_job.mock_calls[count][1][0]
#             self.assertEqual(job.id, job_called_at_count.id)

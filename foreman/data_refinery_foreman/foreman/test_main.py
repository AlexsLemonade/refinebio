from unittest.mock import patch, MagicMock
from datetime import timedelta
from django.utils import timezone
from django.test import TestCase
from data_refinery_foreman.foreman import main
from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    ProcessorJob,
    OriginalFile,
    Dataset,
    SurveyJobKeyValue,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation,
    ProcessorJobDatasetAssociation,
)
from test.support import EnvironmentVarGuard # Python >=3

class ForemanTestCase(TestCase):
    def create_downloader_job(self):
        job = DownloaderJob(downloader_task="SRA",
                            nomad_job_id="DOWNLOADER/dispatch-1528945054-e8eaf540",
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_failed_downloader_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_downloader_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_downloader_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_downloader_jobs.__wrapped__()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = DownloaderJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_lost_downloader_jobs_time(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_downloader_job()
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + timedelta(minutes=1))
        job.save()

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_downloader_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_downloader_jobs.__wrapped__()
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

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_processor_job(self, mock_send_job):
        mock_send_job.return_value = True

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

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_requeuing_processor_job_w_more_ram(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job(pipeline="SALMON", ram_amount=8192)

        main.requeue_processor_job(job)
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)
        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)
        self.assertEqual(original_job.ram_amount, 8192)
        self.assertEqual(retried_job.ram_amount, 12288)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_repeated_processor_failures(self, mock_send_job):
        mock_send_job.return_value = True

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

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_failed_processor_jobs(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job()
        job.success = False
        job.save()

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_failed_processor_jobs.__wrapped__()
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
    def test_retrying_hung_processor_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_processor_jobs.__wrapped__()
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
    def test_not_retrying_hung_processor_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_processor_jobs.__wrapped__()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    @patch('data_refinery_foreman.foreman.main.Nomad')
    def test_retrying_lost_processor_jobs(self, mock_nomad, mock_send_job):
        mock_send_job.return_value = True

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_processor_jobs.__wrapped__()

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
    def test_not_retrying_lost_processor_jobs(self, mock_nomad, mock_send_job):
        """Make sure that we don't retry processor jobs we shouldn't."""
        mock_send_job.return_value = True

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_processor_jobs.__wrapped__()
        self.assertEqual(len(mock_send_job.mock_calls), 0)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertFalse(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertEqual(original_job.success, None)

        # Make sure no additional job was created.
        self.assertEqual(jobs.count(), 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_retrying_lost_processor_jobs_time(self, mock_send_job):
        mock_send_job.return_value = True

        job = self.create_processor_job()
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + timedelta(minutes=1))
        job.save()

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_processor_jobs.__wrapped__()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = ProcessorJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_failed_survey_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_survey_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_hung_survey_jobs.__wrapped__()
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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_survey_jobs.__wrapped__()

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

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_survey_jobs.__wrapped__()
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
        job.created_at = timezone.now() - (main.MIN_LOOP_TIME + timedelta(minutes=1))
        job.save()

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.retry_lost_survey_jobs.__wrapped__()
        self.assertEqual(len(mock_send_job.mock_calls), 1)

        jobs = SurveyJob.objects.order_by('id')
        original_job = jobs[0]
        self.assertTrue(original_job.retried)
        self.assertEqual(original_job.num_retries, 0)
        self.assertFalse(original_job.success)

        retried_job = jobs[1]
        self.assertEqual(retried_job.num_retries, 1)

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_max_jobs(self, mock_send_job):
        return

    @patch('data_refinery_foreman.foreman.main.send_job')
    def test_janitor(self, mock_send_job):
        mock_send_job.return_value = True

        for p in ["1", "2", "3"]:
            pj = ProcessorJob()
            pj.volume_index = p
            pj.save()

        # Just run it once, not forever so get the function that is
        # decorated with @do_forever
        main.send_janitor_jobs.__wrapped__()

        self.assertEqual(ProcessorJob.objects.all().count(), 6)
        self.assertEqual(ProcessorJob.objects.filter(pipeline_applied="JANITOR").count(), 3)

        # Make sure that the janitors are dispatched to the correct volumes.
        ixs = ["1", "2", "3"]
        for p in ProcessorJob.objects.filter(pipeline_applied="JANITOR"):
            self.assertTrue(p.volume_index in ixs)
            ixs.remove(p.volume_index)

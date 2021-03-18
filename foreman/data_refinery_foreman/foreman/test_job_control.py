import datetime
from unittest.mock import patch

from django.test import TestCase, TransactionTestCase

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
)
from data_refinery_foreman.foreman import job_control, utils

# For use in tests that test the JOB_CREATED_AT_CUTOFF functionality.
DAY_BEFORE_JOB_CUTOFF = utils.JOB_CREATED_AT_CUTOFF - datetime.timedelta(days=1)

EMPTY_JOB_QUEUE_RESPONSE = {"jobSummaryList": []}


class ForemanTestCase(TestCase):
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

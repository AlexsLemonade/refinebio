from unittest.mock import patch
from django.test import TestCase
import datetime
from data_refinery_common.models import SurveyJob
from data_refinery_foreman.surveyor import surveyor


class RunJobTestCase(TestCase):
    def test_run_unsupported_source(self):
        """If source_type is unsupported the job still is started and ended."""
        job = SurveyJob(source_type="UNSUPPORTED")
        job.save()
        surveyor.run_job(job)

        self.assertIsInstance(job.replication_ended_at, datetime.datetime)
        self.assertIsInstance(job.start_time, datetime.datetime)
        self.assertIsInstance(job.end_time, datetime.datetime)
        self.assertFalse(job.success)

    @patch.object(surveyor.ArrayExpressSurveyor, "survey")
    def test_calls_survey(self, survey_method):
        """If source_type is supported calls the appropriate survey method."""
        survey_method.return_value = True

        job = SurveyJob(source_type="ARRAY_EXPRESS")
        job.save()
        surveyor.run_job(job)

        survey_method.assert_called()
        self.assertIsInstance(job.replication_ended_at, datetime.datetime)
        self.assertIsInstance(job.start_time, datetime.datetime)
        self.assertIsInstance(job.end_time, datetime.datetime)
        self.assertTrue(job.success)

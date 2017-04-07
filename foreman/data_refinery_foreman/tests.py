from django.test import TestCase
import datetime
from django.utils import timezone
from data_refinery_models.models import SurveyJob
from data_refinery_foreman.surveyor import surveyor


class RunJobTestCase(TestCase):
    def test_run_unsupported_source(self):
        """If source_type is unsupported the job still is started and ended"""
        job = SurveyJob(source_type="UNSUPPORTED")
        surveyor.run_job(job)
        # type(job.replication_ended_at)
        self.assertIsInstance(job.replication_ended_at, datetime.datetime)

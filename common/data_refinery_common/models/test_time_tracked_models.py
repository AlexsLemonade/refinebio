from django.test import TestCase
from django.utils import timezone
from data_refinery_common.models import SurveyJob


class TimeTrackedModelTestCase(TestCase):
    def test_time_tracking_works(self):
        """created_at and updated_at are initialized upon creation"""
        job = SurveyJob.objects.create(source_type="ARRAY_EXPRESS")
        timedelta = timezone.now() - job.created_at
        self.assertLess(timedelta.total_seconds(), 1)
        self.assertEqual(job.created_at, job.updated_at)

        # When the job is updated and saved, updated_at changes
        job.success = False
        job.save()
        self.assertNotEqual(job.created_at, job.updated_at)

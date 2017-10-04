from django.test import TestCase
from django.utils import timezone
from data_refinery_common.models import SurveyJob, Batch, File


class TimeTrackedModelTestCase(TestCase):
    def setUp(self):
        # Fix for issue explained here:
        # https://stackoverflow.com/questions/31504591/interfaceerror-connection-already-closed-using-django-celery-scrapy # noqa
        File.objects.all().delete()
        Batch.objects.all().delete()
        SurveyJob.objects.all().delete()
        # Using this hackish fix because I'll be ripping celery out of
        # the project at some point so a more thorough fix would be a
        # waste of time

        SurveyJob.objects.create(source_type="ARRAY_EXPRESS")

    def tearDown(self):
        SurveyJob.objects.all().delete()

    def test_time_tracking_works(self):
        """created_at and updated_at are initialized upon creation"""
        job = SurveyJob.objects.get(source_type="ARRAY_EXPRESS")
        timedelta = timezone.now() - job.created_at
        self.assertLess(timedelta.total_seconds(), 1)
        self.assertEqual(job.created_at, job.updated_at)

        # When the job is updated and saved, updated_at changes
        job.success = False
        job.save()
        self.assertNotEqual(job.created_at, job.updated_at)

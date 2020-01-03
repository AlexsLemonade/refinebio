from django.test import TestCase
from django.utils import timezone

from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    ProcessorJob,
)


class SanityTestJobsTestCase(TestCase):
    def test_jobs_sanity(self):
        """Just makes sure creating Jobs doesn't fail"""

        s_job = SurveyJob()
        s_job.save()

        processor_job = ProcessorJob()
        processor_job.pipeline_applied = "test0"
        processor_job.save()

        dl_job = DownloaderJob()
        dl_job.downloader_task = "XYZ"
        dl_job.accession_code = "123"
        dl_job.save()

    def test_time_tracking_works(self):
        """created_at and last_modified are initialized upon creation"""
        job = SurveyJob.objects.create(source_type="ARRAY_EXPRESS")
        timedelta = timezone.now() - job.created_at
        self.assertLess(timedelta.total_seconds(), 1)
        self.assertEqual(job.created_at, job.last_modified)

        # When the job is updated and saved, last_modified changes
        job.success = False
        job.save()
        self.assertNotEqual(job.created_at, job.last_modified)

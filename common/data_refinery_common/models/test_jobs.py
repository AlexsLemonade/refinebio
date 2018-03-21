from django.test import TestCase
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
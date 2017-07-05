from unittest.mock import Mock, patch
from django.test import TestCase
from data_refinery_foreman.foreman import main
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    DownloaderJob,
    ProcessorJob
)


class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    def tearDown(self):
        Batch.objects.all().delete()
        DownloaderJob.objects.all().delete()
        ProcessorJob.objects.all().delete()

    def insert_batch(self):
        download_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        batch = Batch(
            survey_job=self.survey_job,
            source_type="ARRAY_EXPRESS",
            size_in_bytes=0,
            download_url=download_url,
            raw_format="CEL",
            processed_format="PCL",
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-1",
            experiment_accession_code="E-MTAB-3050",
            experiment_title="It doesn't really matter.",
            name="CE1234.CEL",
            internal_location="A-AFFY-1/AFFY_TO_PCL/",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch.save()
        return batch

    @patch('data_refinery_foreman.surveyor.message_queue.app.send_task')
    def test_requeuing_downloader_job(self, mock_send_task):
        batch = self.insert_batch()
        job = DownloaderJob.create_job_and_relationships(
            num_retries=0, batches=[batch], downloader_task="dummy_task")

        main.requeue_downloader_job(job)
        mock_send_task.assert_called_once()

        jobs = DownloaderJob.objects.order_by('id')
        self.assertTrue(jobs[0].retried)
        self.assertEqual(jobs[0].num_retries, 0)
        self.assertFalse(jobs[0].success)

        self.assertEqual(jobs[1].num_retries, 1)

import copy
from unittest.mock import patch, MagicMock
from django.test import TestCase
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_workers.downloaders import array_express


class DownloadArrayExpressTestCase(TestCase):
    def test_good_batch_grouping(self):
        """Returns true if all batches have the same download_url."""
        batches = [Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com")]
        job_id = 1
        self.assertIsNone(array_express._verify_batch_grouping(batches, job_id))

    def test_bad_batch_grouping(self):
        """Raises exception if all batches don't have the same download_url."""
        batches = [Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com"),
                   Batch(download_url="https://wompwomp.com")]
        job_id = 1
        with self.assertRaises(ValueError):
            array_express._verify_batch_grouping(batches, job_id)

    @patch("data_refinery_workers.downloaders.array_express.utils.processor_pipeline_registry")
    @patch("data_refinery_workers.downloaders.array_express._verify_batch_grouping")
    @patch("data_refinery_workers.downloaders.array_express._download_file")
    @patch("data_refinery_workers.downloaders.array_express._extract_file")
    def test_download(self,
                      _extract_file,
                      _download_file,
                      _verify_batch_grouping,
                      pipeline_registry):
        # Set up mocks:
        mock_processor_task = MagicMock()
        mock_processor_task.delay = MagicMock()
        mock_processor_task.delay.return_value = None
        pipeline_registry.__getitem__ = MagicMock()
        pipeline_registry.__getitem__.return_value = mock_processor_task

        # Set up database records:
        survey_job = SurveyJob(
            source_type="ARRAY_EXPRESS"
        )
        survey_job.save()

        download_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        batch = Batch(
            survey_job=survey_job,
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
        batch2 = copy.deepcopy(batch)
        batch2.name = "CE2345.CEL"
        batch.save()
        batch2.save()

        downloader_job = DownloaderJob.create_job_and_relationships(batches=[batch, batch2])

        # Call the task we're testing:
        array_express.download_array_express.apply(args=(downloader_job.id,)).get()

        target_file_path = ("/home/user/data_store/temp/A-AFFY-1/AFFY_TO_PCL/{}"
                            "/E-GEOD-59071.raw.3.zip").format(str(downloader_job.id))

        # Verify that all expected functionality is run:
        _verify_batch_grouping.assert_called_once()
        _download_file.assert_called_with(download_url, target_file_path, downloader_job.id)
        _extract_file.assert_called_with([batch, batch2], downloader_job.id)

        mock_processor_task.delay.assert_called()

        # Verify that the database has been updated correctly:
        batches = Batch.objects.all()
        for batch in batches:
            self.assertEqual(batch.status, BatchStatuses.DOWNLOADED.value)

        downloader_job = DownloaderJob.objects.get()
        self.assertTrue(downloader_job.success)
        self.assertIsNotNone(downloader_job.end_time)

        processor_jobs = ProcessorJob.objects.all()
        self.assertEqual(len(processor_jobs), 2)

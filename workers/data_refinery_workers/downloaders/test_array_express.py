import copy
from typing import List
from unittest.mock import patch, call
from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_workers.downloaders import array_express, utils
from data_refinery_common.job_lookup import ProcessorPipeline


class DownloadArrayExpressTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    def insert_objects(self) -> List[Batch]:
        download_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        batch = Batch(
            survey_job=self.survey_job,
            source_type="ARRAY_EXPRESS",
            pipeline_required="AFFY_TO_PCL",
            platform_accession_code="A-AFFY-1",
            experiment_accession_code="E-MTAB-3050",
            experiment_title="It doesn't really matter.",
            organism_id=9606,
            organism_name="HOMO SAPIENS",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch2 = copy.deepcopy(batch)
        batch.save()
        batch2.save()

        file = File(size_in_bytes=0,
                    download_url=download_url,
                    raw_format="CEL",
                    processed_format="PCL",
                    name="CE1234.CEL",
                    internal_location="A-AFFY-1/AFFY_TO_PCL/",
                    batch=batch)
        file2 = File(size_in_bytes=0,
                     download_url=download_url,
                     raw_format="CEL",
                     processed_format="PCL",
                     name="CE2345.CEL",
                     internal_location="A-AFFY-1/AFFY_TO_PCL/",
                     batch=batch2)
        file.save()
        file2.save()

        batch.files = [file]
        batch2.files = [file]

        return ([batch, batch2], [file, file2])

    def test_good_batch_grouping(self):
        """Returns true if all batches have the same download_url."""
        batches, files = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(
            batches=batches, downloader_task="dummy")

        self.assertIsNone(array_express._verify_batch_grouping(files, downloader_job))

    def test_bad_batch_grouping(self):
        """Raises exception if all batches don't have the same download_url."""
        batches, files = self.insert_objects()
        files[1].download_url = "https://wompwomp.com"
        files[1].save()
        downloader_job = DownloaderJob.create_job_and_relationships(
            batches=batches, downloader_task="dummy")

        with self.assertRaises(ValueError):
            array_express._verify_batch_grouping(files, downloader_job)

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch("data_refinery_workers.downloaders.array_express._verify_batch_grouping")
    @patch("data_refinery_workers.downloaders.array_express._download_file")
    @patch("data_refinery_workers.downloaders.array_express._extract_file")
    def test_download(self,
                      _extract_file,
                      _download_file,
                      _verify_batch_grouping,
                      mock_send_job):
        mock_send_job.return_value = None

        batches, files = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the function we're testing:
        array_express.download_array_express(downloader_job.id)

        target_file_path = ("/home/user/data_store/temp/A-AFFY-1/AFFY_TO_PCL/downloader_job_{}"
                            "/E-GEOD-59071.raw.3.zip").format(str(downloader_job.id))

        # Verify that all expected functionality is run:
        self.assertEqual(_verify_batch_grouping.call_count, 1)
        download_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        _download_file.assert_called_with(download_url, target_file_path, downloader_job)
        args, _ = _extract_file.call_args
        file_query_set, job = args
        self.assertEqual(list(file_query_set), files)
        self.assertEqual(job.id, downloader_job.id)

        # Verify that the database has been updated correctly:
        batches = Batch.objects.all()
        for batch in batches:
            self.assertEqual(batch.status, BatchStatuses.DOWNLOADED.value)

        downloader_job = DownloaderJob.objects.get()
        self.assertTrue(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)

        processor_jobs = ProcessorJob.objects.all()
        self.assertEqual(len(processor_jobs), 2)

        mock_send_job.assert_has_calls([
            call(ProcessorPipeline.AFFY_TO_PCL,
                 processor_jobs[0].id),
            call(ProcessorPipeline.AFFY_TO_PCL,
                 processor_jobs[1].id)
        ])

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch("data_refinery_workers.downloaders.array_express._download_file")
    @patch("data_refinery_workers.downloaders.array_express._extract_file")
    def test_verification_failure(self,
                                  _extract_file,
                                  _download_file,
                                  mock_send_job):
        mock_send_job.return_value = None

        # Set a different download URL to trigger a failure in the
        # _verify_batch_grouping function
        batches, files = self.insert_objects()
        files[1].download_url = "https://wompwomp.com"
        files[1].save()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function
        array_express.download_array_express(downloader_job.id)

        _download_file.assert_not_called()
        _extract_file.assert_not_called()
        mock_send_job.assert_not_called()

        # Verify that the database has been updated correctly:
        downloader_job = DownloaderJob.objects.get()
        self.assertFalse(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)
        self.assertEqual(downloader_job.failure_reason,
                         ("A Batch's file doesn't have the same download "
                          "URL as the other batches' files."))

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch('data_refinery_workers.downloaders.array_express.open')
    @patch("data_refinery_workers.downloaders.array_express._extract_file")
    def test_download_failure(self,
                              _extract_file,
                              _open,
                              mock_send_job):
        # Set up mocks:
        mock_send_job.return_value = None
        _open.side_effect = Exception()

        batches, _ = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function
        array_express.download_array_express(downloader_job.id)

        _extract_file.assert_not_called()
        mock_send_job.assert_not_called()

        # Verify that the database has been updated correctly:
        downloader_job = DownloaderJob.objects.get()
        self.assertFalse(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)
        self.assertEqual(downloader_job.failure_reason,
                         "Exception caught while downloading batch")

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch("data_refinery_workers.downloaders.array_express._download_file")
    def test_extraction_failure(self,
                                _download_file,
                                mock_send_job):
        mock_send_job.return_value = None

        batches, files = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function
        array_express.download_array_express(downloader_job.id)

        mock_send_job.assert_not_called()

        # Verify that the database has been updated correctly:
        downloader_job = DownloaderJob.objects.get()
        self.assertFalse(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)

        job_dir = utils.JOB_DIR_PREFIX + str(downloader_job.id)
        zip_path = files[0].get_temp_download_path(job_dir)
        self.assertEqual(downloader_job.failure_reason,
                         "Exception caught while extracting " + zip_path)

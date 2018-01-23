import shutil
import copy
from typing import List
from unittest.mock import patch, call
from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchKeyValue,
    BatchStatuses,
    File,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_workers.downloaders import transcriptome_index
from data_refinery_common.job_lookup import ProcessorPipeline


class DownloadTranscriptomeIndexTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job
        self.gtf_download_url = ("ftp://ftp.ensemblgenomes.org/pub/release-37/plants/gtf/"
                                 "aegilops_tauschii/Aegilops_tauschii.ASM34733v1.37.gtf.gz")
        self.fasta_download_url = (
            "ftp://ftp.ensemblgenomes.org/pub/release-37/plants/fasta/"
            "aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz")

    def insert_objects(self) -> List[Batch]:
        download_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"  # noqa
        batch1 = Batch(
            survey_job=self.survey_job,
            source_type="TRANSCRIPTOME_INDEX",
            pipeline_required="TRANSCRIPTOME_INDEX",
            platform_accession_code="EnsemblPlants",
            experiment_accession_code="AEGILOPS_TAUSCHII",
            experiment_title="It doesn't really matter.",
            organism_id=37682,
            organism_name="AEGILOPS TAUSCHII",
            release_date="2017-05-05",
            last_uploaded_date="2017-05-05",
            status=BatchStatuses.NEW.value
        )
        batch2 = copy.deepcopy(batch1)
        batch1.save()
        batch2.save()

        for batch, length, kmer_size in [(batch1, "_short", "23"), (batch2, "_long", "31")]:
            BatchKeyValue(batch=batch, key="length", value=length).save()
            BatchKeyValue(batch=batch, key="kmer_size", value=kmer_size).save()

            file1 = File(size_in_bytes=0,
                         download_url=self.fasta_download_url,
                         raw_format="fa.gz",
                         processed_format="tar.gz",
                         name="Aegilops_tauschii{}.fa.gz".format(length),
                         internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
                         batch=batch)
            file2 = File(size_in_bytes=0,
                         download_url=self.gtf_download_url,
                         raw_format="gtf.gz",
                         processed_format="tar.gz",
                         name="Aegilops_tauschii{}.gtf.gz".format(length),
                         internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
                         batch=batch)
            file1.save()
            file2.save()
            batch.files = [file1, file2]

        return [batch1, batch2]

    def test_good_file_grouping(self):
        """Returns None if both files have the same download_url."""
        downloader_job = DownloaderJob.create_job_and_relationships(
            batches=[], downloader_task="dummy")

        self.assertIsNone(transcriptome_index._verify_files(File(download_url="a"),
                                                            File(download_url="a"),
                                                            downloader_job))

    def test_bad_file_grouping(self):
        """Raises exception if both files don't have the same download_url."""
        downloader_job = DownloaderJob.create_job_and_relationships(
            batches=[], downloader_task="dummy")

        with self.assertRaises(ValueError):
            self.assertIsNone(transcriptome_index._verify_files(File(download_url="a"),
                                                                File(download_url="b"),
                                                                downloader_job))

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch("data_refinery_workers.downloaders.transcriptome_index._verify_files")
    @patch("data_refinery_workers.downloaders.transcriptome_index._download_file")
    @patch("data_refinery_workers.downloaders.transcriptome_index._upload_files")
    def test_download(self,
                      _upload_files,
                      _download_file,
                      _verify_files,
                      mock_send_job):
        # Clean up temp directory:
        shutil.rmtree("/home/user/data_store/temp/EnsemblPlants/TRANSCRIPTOME_INDEX",
                      ignore_errors=True)

        mock_send_job.return_value = None

        batches = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function we're testing:
        transcriptome_index.download_transcriptome(downloader_job.id)

        target_gtf_path = (
            "/home/user/data_store/temp/EnsemblPlants/TRANSCRIPTOME_INDEX/downloader_job_{}"
            "/Aegilops_tauschii_short.gtf.gz"
        ).format(str(downloader_job.id))
        target_fasta_path = (
            "/home/user/data_store/temp/EnsemblPlants/TRANSCRIPTOME_INDEX/downloader_job_{}"
            "/Aegilops_tauschii_short.fa.gz"
        ).format(str(downloader_job.id))

        # Verify that all expected functionality is run:
        self.assertEqual(_verify_files.call_count, 2)
        self.assertEqual(_download_file.call_count, 2)
        _download_file.assert_any_call(self.gtf_download_url, target_gtf_path, downloader_job)
        _download_file.assert_any_call(self.fasta_download_url, target_fasta_path, downloader_job)
        args, _ = _upload_files.call_args
        job_dir, files, job = args
        self.assertEqual(set(files), set(batches[0].files + batches[1].files))
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
            call(ProcessorPipeline.TRANSCRIPTOME_INDEX,
                 processor_jobs[0].id),
            call(ProcessorPipeline.TRANSCRIPTOME_INDEX,
                 processor_jobs[1].id)
        ])

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch("data_refinery_workers.downloaders.transcriptome_index._download_file")
    @patch("data_refinery_workers.downloaders.transcriptome_index._upload_files")
    def test_verification_failure(self, _upload_files, _download_file, mock_send_job):
        mock_send_job.return_value = None

        # Set a different download URL to trigger a failure in the
        # _verify_batch_grouping function
        batches = self.insert_objects()
        batches[0].files[0].download_url = "https://wompwomp.com"
        batches[0].files[0].save()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function
        transcriptome_index.download_transcriptome(downloader_job.id)

        _download_file.assert_not_called()
        _upload_files.assert_not_called()
        mock_send_job.assert_not_called()

        # Verify that the database has been updated correctly:
        downloader_job = DownloaderJob.objects.get()
        self.assertFalse(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)
        self.assertEqual(downloader_job.failure_reason,
                         ("A Batch's file doesn't have the same download "
                          "URL as the other batch's file."))

    @patch("data_refinery_workers.downloaders.utils.send_job")
    @patch('data_refinery_workers.downloaders.transcriptome_index.open')
    @patch("data_refinery_workers.downloaders.transcriptome_index._upload_files")
    def test_download_failure(self,
                              _upload_files,
                              _open,
                              mock_send_job):
        # Set up mocks:
        mock_send_job.return_value = None
        _open.side_effect = Exception()

        batches = self.insert_objects()
        downloader_job = DownloaderJob.create_job_and_relationships(batches=batches)

        # Call the downloader function
        transcriptome_index.download_transcriptome(downloader_job.id)

        _upload_files.assert_not_called()
        mock_send_job.assert_not_called()

        # Verify that the database has been updated correctly:
        downloader_job = DownloaderJob.objects.get()
        self.assertFalse(downloader_job.success)
        self.assertIsNotNone(downloader_job.start_time)
        self.assertIsNotNone(downloader_job.end_time)
        failure_reason = "Exception caught while downloading file from: {}".format(
            batches[0].files[0].download_url)
        self.assertEqual(downloader_job.failure_reason, failure_reason)

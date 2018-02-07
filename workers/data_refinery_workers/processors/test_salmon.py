import os
import shutil
from django.test import TestCase, tag
from unittest.mock import patch
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchKeyValue,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import salmon, utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


def init_objects():
    survey_job = SurveyJob(source_type="SALMON")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="SALMON",
        pipeline_required="SALMON",
        platform_accession_code="IlluminaGenomeAnalyzerII",
        experiment_accession_code="ERX000259",
        experiment_title="It doesn't really matter.",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-11-02",
        last_uploaded_date="2017-11-02",
        status=BatchStatuses.DOWNLOADED.value
    )
    batch.save()

    first_fastq_file = File(
        size_in_bytes=2214725074,
        raw_format="fastq.gz",
        processed_format="tar.gz",
        name="ERR003000_1.fastq.gz",
        internal_location="IlluminaGenomeAnalyzerII/SALMON",
        download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR003/"
                      "ERR003000/ERR003000_1.fastq.gz"),
        batch=batch
    )
    first_fastq_file.save()

    second_fastq_file = File(
        size_in_bytes=2214725074,
        raw_format="fastq.gz",
        processed_format="tar.gz",
        name="ERR003000_2.fastq.gz",
        internal_location="IlluminaGenomeAnalyzerII/SALMON",
        download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR003/"
                      "ERR003000/ERR003000_2.fastq.gz"),
        batch=batch
    )
    second_fastq_file.save()

    batch.files = [first_fastq_file, second_fastq_file]
    return (batch, first_fastq_file, second_fastq_file)


def _insert_salmon_index():
    """Creates a batch for the index for the organism for the test."""
    survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="TRANSCRIPTOME_INDEX",
        pipeline_required="TRANSCRIPTOME_INDEX",
        platform_accession_code="TEST",
        experiment_accession_code="HOMO_SAPIENS",
        experiment_title="It doesn't really matter.",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-11-02",
        last_uploaded_date="2017-11-02",
        status=BatchStatuses.PROCESSED.value
    )
    batch.save()

    kmer_size = BatchKeyValue(
        key="kmer_size",
        value="23",
        batch=batch
    )
    kmer_size.save()

    index_file = File(
        size_in_bytes=2214725074,
        raw_format="gtf.gz",
        processed_format="tar.gz",
        name="Homo_sapiens_short.gtf.gz",
        internal_location="TEST/TRANSCRIPTOME_INDEX",
        download_url=("ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens"
                      "/Homo_sapiens.GRCh38.90.gtf.gz"),
        batch=batch
    )
    index_file.save()


class SalmonTestCase(TestCase):
    """Test all functionality for the SALMON processor.
    The successful functionality of all individual functions are

    tested in test_success since each function sets up the next one
    nicely. The failures of all individual functions are tested in
    separate test functions.
    """

    @tag('slow')
    def test_success(self):
        """Tests the successful path of the module under test."""
        logger.info("STARTING SALMON SUCCESS TEST!!!!!!!!")
        # Set up test environment.
        batch, first_file, second_file = init_objects()
        _insert_salmon_index()
        # Change the batch/files to point to test-specific locations
        batch.platform_accession_code = "TEST"
        batch.save()
        first_file.internal_location = "TEST/SALMON"
        first_file.save()
        second_file.internal_location = "TEST/SALMON"
        second_file.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        processor_job.save()
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = salmon._set_job_prefix(job_context)

        # Ensure temp dir isn't leftover from a previous test.
        temp_dir = first_file.get_temp_dir(job_context["job_dir_prefix"])
        shutil.rmtree(temp_dir, ignore_errors=True)

        # One of the functions being tested:
        job_context = salmon._prepare_files(job_context)

        input_file_path = job_context["input_file_path"]
        self.assertIsInstance(input_file_path, str)
        self.assertTrue(os.path.isfile(input_file_path))
        input_file_path_2 = job_context["input_file_path_2"]
        self.assertIsInstance(input_file_path_2, str)
        self.assertTrue(os.path.isfile(input_file_path_2))
        output_directory_path = job_context["output_directory"]
        self.assertIsInstance(output_directory_path, str)
        self.assertTrue(os.path.isdir(output_directory_path))

        job_context = salmon._determine_index_length(job_context)

        # The 'kmer_size' key has been added to job_context with the
        # correct value.
        self.assertEqual(job_context["kmer_size"], "23")

        # Another function being tested
        job_context = salmon._download_index(job_context)

        self.assertTrue(job_context["success"])
        self.assertTrue("index_directory" in job_context)
        self.assertTrue(os.path.isdir(job_context["index_directory"]))
        self.assertEqual(9, len(os.listdir(job_context["index_directory"])))

        # Another function being tested
        job_context = salmon._run_salmon(job_context)

        self.assertTrue(job_context["success"])
        self.assertGreater(len(os.listdir(output_directory_path)), 1)

        # The last function to test
        job_context = salmon._zip_and_upload(job_context)

        self.assertTrue(job_context["success"])
        self.assertTrue(os.path.exists(first_file.get_processed_path()))

        # Clean up both input and output files
        first_file.remove_temp_directory()
        shutil.rmtree(first_file.get_processed_dir())
        logger.info("ENDING SALMON SUCCESS TEST!!!!!!!!")

    def test_prepare_files_failure(self):
        batch, _, _ = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = salmon._prepare_files(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while retrieving raw file ERR003000_2.fastq.gz")

        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

    def test_download_index_not_found(self):
        batch, _, _ = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "kmer_size": "23"})
        job_context = salmon._prepare_files(job_context)

        # Function we're testing.
        salmon._download_index(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Failed to find an index for organism HOMO SAPIENS with kmer_size of 23.")

    @patch.object(File, "download_processed_file")
    def test_download_index_missing(self, mock_download_processed_file):
        batch, _, _ = init_objects()
        _insert_salmon_index()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "job_dir_prefix": "dummy",
                                       "kmer_size": "23"})
        job_context = salmon._prepare_files(job_context)

        mock_download_processed_file.side_effect = FileNotFoundError()

        # The function being testing.
        salmon._download_index(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Failed to download and extract index tarball Homo_sapiens_short.gtf.gz")

    def test_run_salmon_failure(self):
        batch, first_fastq_file, second_fastq_file = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # Mock out the job_context with everything the function under
        # test will expect
        input_file_path_1 = first_fastq_file.get_temp_pre_path("dummy")
        input_file_path_2 = second_fastq_file.get_temp_pre_path("dummy")
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "job_dir_prefix": "dummy",
                                       "batches": [batch],
                                       "index_directory": "missing",
                                       "input_file_path": input_file_path_1,
                                       "input_file_path_2": input_file_path_2,
                                       "output_directory": "blah"})

        # The function being tested.
        job_context = salmon._run_salmon(job_context)

        self.assertFalse(job_context["success"])
        self.assertNotEqual(processor_job.failure_reason, None)
        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

    def test_zip_and_upload_failure(self):
        # Initialize test objects
        batch, _, _ = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # Mock out the job_context with everything the function under
        # test will expect
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "job_dir_prefix": "dummy",
                                       "output_directory": "missing/index"})

        # The function being tested.
        job_context = salmon._zip_and_upload(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while zipping processed directory ERR003000_1.fastq.gz")
        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

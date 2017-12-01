import os
import shutil
from django.test import TestCase
from unittest.mock import patch
from subprocess import CompletedProcess
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchKeyValue,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import salmon, utils


def init_objects():
    survey_job = SurveyJob(source_type="SALMON")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="SALMON",
        pipeline_required="SALMON",
        platform_accession_code="IlluminaHiSeq2500",
        experiment_accession_code="PRJEB5019",
        experiment_title="It doesn't really matter.",
        organism_id=100900,
        organism_name="MUS MUSCULUS",
        release_date="2017-11-02",
        last_uploaded_date="2017-11-02",
        status=BatchStatuses.DOWNLOADED.value
    )
    batch.save()

    first_fastq_file = File(
        size_in_bytes=2214725074,
        raw_format="fastq",
        processed_format="sf",
        name="ERR1680082_1.fastq",
        internal_location="IlluminaHiSeq2500/SALMON",
        download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/"
                      "ERR1680082/ERR1680082_1.fastq.gz"),
        batch=batch
    )
    first_fastq_file.save()

    second_fastq_file = File(
        size_in_bytes=2214725074,
        raw_format="fastq",
        processed_format="sf",
        name="ERR1680082_2.fastq",
        internal_location="IlluminaHiSeq2500/SALMON",
        download_url=("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168/002/"
                      "ERR1680082/ERR1680082_2.fastq.gz"),
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
        experiment_accession_code="MUS_MUSCULUS",
        experiment_title="It doesn't really matter.",
        organism_id=100900,
        organism_name="MUS MUSCULUS",
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
        name="Mus_musculus_short.gtf.gz",
        internal_location="TEST/TRANSCRIPTOME_INDEX",
        download_url=("ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus"
                      "/Mus_musculus.GRCm38.90.gtf.gz"),
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

    def test_success(self):
        """Tests the successful path of the module under test."""
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

        self.assertFalse("success" in job_context)
        self.assertTrue("index_directory" in job_context)
        self.assertTrue(os.path.isdir(job_context["index_directory"]))
        self.assertEqual(9, len(os.listdir(job_context["index_directory"])))

        # Another function being tested
        job_context = salmon._run_salmon(job_context)

        self.assertFalse("success" in job_context)
        self.assertEqual(6, len(os.listdir(output_directory_path)))

        # The last function to test
        job_context = salmon._zip_and_upload(job_context)

        self.assertTrue(job_context["success"])
        self.assertTrue(os.path.exists(first_file.get_processed_path()))

        # Clean up both input and output files
        first_file.remove_temp_directory()
        shutil.rmtree(first_file.get_processed_dir())

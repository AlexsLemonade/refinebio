import os
import shutil
from django.test import TestCase
from unittest.mock import patch
from subprocess import CompletedProcess
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import transcriptome_index, utils


def init_objects():
    survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="TRANSCRIPTOME_INDEX",
        pipeline_required="TRANSCRIPTOME_INDEX",
        platform_accession_code="EnsemblPlants",
        experiment_accession_code="aegilops_tauschii",
        experiment_title="It doesn't really matter.",
        organism_id=37682,
        organism_name="AEGILOPS TAUSCHII",
        release_date="2017-11-02",
        last_uploaded_date="2017-11-02",
        status=BatchStatuses.DOWNLOADED.value
    )
    batch.save()

    gtf_file = File(
        size_in_bytes=-1,
        raw_format="gtf.gz",
        processed_format="tar.gz",
        name="aegilops_tauschii_short.gtf.gz",
        internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
        download_url=("ftp://ftp.ensemblgenomes.org/pub/release-37/plants/gtf"
                      "/aegilops_tauschii/Aegilops_tauschii.ASM34733v1.37.gtf.gz"),
        batch=batch
    )
    gtf_file.save()

    fasta_file = File(
        size_in_bytes=-1,
        raw_format="fa.gz",
        processed_format="tar.gz",
        name="aegilops_tauschii_short.fa.gz",
        internal_location="EnsemblPlants/TRANSCRIPTOME_INDEX",
        download_url=("ftp://ftp.ensemblgenomes.org/pub/release-37/plants/fasta"
                      "/aegilops_tauschii/dna/Aegilops_tauschii."
                      "ASM34733v1.dna.toplevel.fa.gz"),
        batch=batch
    )
    fasta_file.save()

    batch.files = [gtf_file, fasta_file]
    return (batch, gtf_file, fasta_file)


class PrepareFilesTestCase(TestCase):
    """I should switch the testing strategy here from doing each function
    completely seperately to doing them all in the same testing go,
    verifying that each step of the way things are as expected. This
    way I'll have less code duplicating the functionality already
    existing in the file under test.  I should then also test the
    potential failures of each function seperately.
    """

    def test_success(self):
        """Tests the successful path of the module under test."""
        # Set up test environment.
        batch, gtf_file, fasta_file = init_objects()
        # Change the batch/files to point to test-specific locations
        batch.platform_accession_code = "TEST"
        batch.save()
        gtf_file.internal_location = "TEST/TRANSCRIPTOME_INDEX"
        gtf_file.save()
        fasta_file.internal_location = "TEST/TRANSCRIPTOME_INDEX"
        fasta_file.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        processor_job.save()
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = transcriptome_index._set_job_prefix(job_context)

        # Ensure temp dir isn't leftover from a previous test.
        temp_dir = gtf_file.get_temp_dir(job_context["job_dir_prefix"])
        shutil.rmtree(temp_dir, ignore_errors=True)

        # This will be useful in failure tests
        # for file in batch.files:
        #     raw_path = file.get_raw_path()
        #     with gzip.open(raw_path, "wt") as dummy_fastq:
        #         dummy_fastq.write("This is a dummy file for tests to operate upon.")

        # The function being tested:
        job_context = transcriptome_index._prepare_files(job_context)

        gtf_file_path = job_context["gtf_file_path"]
        self.assertIsInstance(gtf_file_path, str)
        self.assertTrue(os.path.isfile(gtf_file_path))
        fasta_file_path = job_context["fasta_file_path"]
        self.assertIsInstance(fasta_file_path, str)
        self.assertTrue(os.path.isfile(fasta_file_path))
        self.assertIsInstance(job_context["gtf_file"], File)
        self.assertIsInstance(job_context["fasta_file"], File)

        job_context = transcriptome_index._process_gtf(job_context)

        # "success" is only populated by this function on an error
        self.assertFalse("success" in job_context)
        # A new gtf file is created by _process_gtf
        self.assertNotEqual(job_context["gtf_file_path"], gtf_file_path)
        self.assertTrue(os.path.isfile(job_context["gtf_file_path"]))
        self.assertTrue("genes_to_transcripts_path" in job_context)
        self.assertTrue(os.path.isfile(job_context["genes_to_transcripts_path"]))

        #
        job_context = transcriptome_index._prepare_reference(job_context)

        self.assertFalse("success" in job_context)
        self.assertTrue("output_dir" in job_context)
        self.assertTrue(os.path.isdir(job_context["output_dir"]))
        # There should be seven output files in the directory
        self.assertEqual(7, len(os.listdir(job_context["output_dir"])))

        job_context = transcriptome_index._zip_index(job_context)

        self.assertTrue("files_to_upload" in job_context)
        self.assertEqual(job_context["gtf_file"].id, gtf_file.id)
        self.assertEqual(job_context["files_to_upload"][0].id, gtf_file.id)
        zipped_path = gtf_file.get_temp_post_path(job_context["job_dir_prefix"])
        self.assertTrue(os.path.isfile(zipped_path))

        # Clean up both input and output files
        job_context["gtf_file"].remove_temp_directory()

    def test_prepare_files_failure(self):
        batch, _, _ = init_objects()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = transcriptome_index._prepare_files(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while retrieving raw files.")

        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

    @patch("data_refinery_workers.processors.transcriptome_index.subprocess.run")
    def test_prepare_reference_failure(self, mocked_subprocess):
        # Initialize mock and test objects
        mocked_subprocess.return_value = CompletedProcess([], 1, stdout=None,
                                                          stderr="Error: something went wrong.")
        batch, gtf_file, fasta_file = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # Mock out the job_context with everything the function under
        # test will expect
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "gtf_file": gtf_file,
                                       "gtf_file_path": "dummy",
                                       "fasta_file": fasta_file,
                                       "fasta_file_path": "dummy",
                                       "genes_to_transcripts_path": "dummy"})
        job_context = transcriptome_index._set_job_prefix(job_context)

        # The function being tested.
        job_context = transcriptome_index._prepare_reference(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         ("Shell call to rsem-prepare-reference failed because: "
                          "Error: something went wrong."))
        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

    def test_zip_index_failure(self):
        # Initialize test objects
        batch, gtf_file, _ = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # Mock out the job_context with everything the function under
        # test will expect
        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id,
                                       "gtf_file": gtf_file,
                                       "output_dir": "missing/index"})
        job_context = transcriptome_index._set_job_prefix(job_context)

        # The function being tested.
        job_context = transcriptome_index._zip_index(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         ("Exception caught while zipping index directory /home/user"
                          "/data_store/temp/EnsemblPlants/TRANSCRIPTOME_INDEX/processor_job_1"
                          "/aegilops_tauschii_short.tar.gz"))
        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))

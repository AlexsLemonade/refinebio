import os
import shutil
from django.test import TestCase
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import salmon, utils


def init_objects():
    survey_job = SurveyJob(source_type="SRA")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="SRA",
        pipeline_required="SALMON",
        platform_accession_code="IlluminaHiSeq2500",
        experiment_accession_code="DRX014494",
        experiment_title="It doesn't really matter.",
        organism_id=10090,
        organism_name="ARABIDOPSIS THALIANA",
        release_date="2015-05-03",
        last_uploaded_date="2015-06-19",
        status=BatchStatuses.DOWNLOADED.value
    )
    batch.save()

    file1 = File(
        size_in_bytes=967794560,
        raw_format="fastq.gz",
        processed_format="tar.gz",
        name="DRR016125_1.fastq.gz",
        internal_location="IlluminaHiSeq2500/SALMON",
        download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_1.fastq.gz",  # noqa
        batch=batch
    )
    file1.save()

    file2 = File(
        size_in_bytes=1001146319,
        raw_format="fastq.gz",
        processed_format="tar.gz",
        name="DRR016125_2.fastq.gz",
        internal_location="IlluminaHiSeq2500/SALMON",
        download_url="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR016125/DRR016125_2.fastq.gz",  # noqa
        batch=batch
    )
    file2.save()

    batch.files = [file1, file2]
    return batch


class PrepareFilesTestCase(TestCase):
    def test_success(self):
        batch = init_objects()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        os.makedirs(batch.files[0].get_raw_dir(), exist_ok=True)
        for file in batch.files:
            raw_path = file.get_raw_path()
            with open(raw_path, "w") as dummy_fastq:
                dummy_fastq.write("This is a dummy file for tests to operate upon.")

        job_context = utils.start_job({"job": processor_job})
        job_context = salmon._prepare_files(job_context)

        input_file_path = job_context["input_file_path"]
        self.assertIsInstance(input_file_path, str)
        self.assertTrue(os.path.isfile(input_file_path))
        input_file_path_2 = job_context["input_file_path_2"]
        self.assertIsInstance(input_file_path_2, str)
        self.assertTrue(os.path.isfile(input_file_path_2))
        self.assertIsInstance(job_context["output_directory"], str)
        os.remove(raw_path)
        os.remove(input_file_path)

    def test_failure(self):
        batch = init_objects()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = salmon._prepare_files(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while retrieving raw file DRR016125_2.fastq.gz")

        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))


class RunSalmonTestCase(TestCase):
    @staticmethod
    def prepare_job_context():
        batch = init_objects()
        # Change the batch/files to point to test-specific locations
        batch.platform_accession_code = "TEST"
        batch.save()
        batch.files[0].internal_location = "TEST/SALMON"
        batch.files[0].save()
        batch.files[1].internal_location = "TEST/SALMON"
        batch.files[1].save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # We have test files in the repo, but they needs to be in the
        # correct location which depends on the ID of the Batch, which
        # changes based on the order tests are run in.
        os.makedirs(batch.files[0].get_temp_dir(), exist_ok=True)
        test_file_path = "/home/user/data_store/temp/TEST/SALMON/" + batch.files[0].name
        input_file_path = batch.files[0].get_temp_pre_path()
        shutil.copyfile(test_file_path, input_file_path)
        test_file_path_2 = "/home/user/data_store/temp/TEST/SALMON/" + batch.files[1].name
        input_file_path_2 = batch.files[1].get_temp_pre_path()
        shutil.copyfile(test_file_path_2, input_file_path_2)

        output_directory = os.path.join(batch.files[0].get_temp_dir(), "output")
        return {"job_id": processor_job.id,
                "job": processor_job,
                "batches": [batch],
                "input_file_path": input_file_path,
                "input_file_path_2": input_file_path_2,
                "index_directory": "/home/user/data_store/mouse_index",
                "output_directory": output_directory}

    def test_success(self):
        job_context = RunSalmonTestCase.prepare_job_context()

        # If output_directory exists, remove it first.
        if os.path.isdir(job_context["output_directory"]):
            shutil.rmtree(job_context["output_directory"])

        job_context = salmon._run_salmon(job_context)

        # success is only populated by this function on an error
        self.assertFalse("success" in job_context)
        self.assertTrue(os.path.isdir(job_context["output_directory"]))

        # Clean up both input and output files
        job_context["batches"][0].files[0].remove_temp_directory()

    def test_failure(self):
        job_context = RunSalmonTestCase.prepare_job_context()
        job_context["index_directory"] = "/home/user/data_store/bad_index"

        job_context = salmon._run_salmon(job_context)
        failure_reason = ("Shell call to salmon failed because: "
                          "Error: The index version file "
                          "/home/user/data_store/bad_index/versionInfo.json doesn't seem to "
                          "exist.  Please try re-building the salmon index.]\\nsalmon quant "
                          "was invoked improperly.\\nFor usage inform")

        self.assertFalse(job_context["success"])
        batch = job_context["batches"][0]
        self.assertFalse(os.path.isfile(batch.files[0].get_temp_dir()))
        self.assertEqual(job_context["job"].failure_reason,
                         failure_reason)

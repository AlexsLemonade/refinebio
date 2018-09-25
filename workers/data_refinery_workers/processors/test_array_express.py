import os
import shutil
from contextlib import closing
from django.test import TestCase
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import array_express, utils
import pandas as pd


def init_objects():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="ARRAY_EXPRESS",
        pipeline_required="AFFY_TO_PCL",
        platform_accession_code="A-AFFY-1",
        experiment_accession_code="E-MTAB-3050",
        experiment_title="It doesn't really matter.",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-05-05",
        last_uploaded_date="2017-05-05",
        status=BatchStatuses.DOWNLOADED.value
    )
    batch.save()

    file = File(size_in_bytes=0,
                download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",  # noqa
                raw_format="CEL",
                processed_format="PCL",
                name="CE1234.CEL",
                internal_location="A-AFFY-1/AFFY_TO_PCL/",
                batch=batch)
    file.save()

    batch.files = [file]
    return batch


class PrepareFilesTestCase(TestCase):
    def test_success(self):
        batch = init_objects()
        file = batch.files[0]

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        os.makedirs(file.get_raw_dir(), exist_ok=True)
        raw_path = file.get_raw_path()
        with open(raw_path, "w") as dummy_pcl:
            dummy_pcl.write("This is a dummy file for tests to operate upon.")

        job_context = utils.start_job({"job": processor_job})
        job_context = array_express._prepare_files(job_context)

        input_file_path = job_context["input_file_path"]
        self.assertIsInstance(input_file_path, str)
        self.assertIsInstance(job_context["output_file_path"], str)

        self.assertTrue(os.path.isfile(input_file_path))
        os.remove(raw_path)
        os.remove(input_file_path)

    def test_failure(self):
        batch = init_objects()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        job_context = utils.start_job({"job": processor_job,
                                       "job_id": processor_job.id})
        job_context = array_express._prepare_files(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while retrieving raw file CE1234.CEL")

        self.assertFalse(os.path.isfile(batch.files[0].get_temp_pre_path()))


class DetermineBrainarrayPackageTestCase(TestCase):
    def test_success(self):
        batch = init_objects()
        file = batch.files[0]
        batch.platform_accession_code = "TEST"
        batch.save()
        file.internal_location = "TEST/AFFY_TO_PCL"
        file.name = "GSM1426186_UC_colon_inactive_201.CEL"
        file.save()

        file.remove_temp_directory = MagicMock()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # We have a test file in the repo, but it needs to be in the
        # correct location which depends on the ID of the Batch, which
        # changes based on the order tests are run in.
        test_file_path = "/home/user/data_store/temp/TEST/AFFY_TO_PCL/" + file.name
        input_file_path = file.get_temp_pre_path()
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        shutil.copyfile(test_file_path, input_file_path)

        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "input_file_path": input_file_path}

        job_context = array_express._determine_brainarray_package(job_context)

        self.assertEqual(job_context["brainarray_package"], "hugene10sthsentrezgprobe")
        file.remove_temp_directory.assert_not_called()

        # Clean up the copied file
        os.remove(input_file_path)

    def test_failure(self):
        batch = init_objects()
        file = batch.files[0]
        batch.platform_accession_code = "TEST2"
        batch.save()
        file.internal_location = "TEST2/AFFY_TO_PCL"
        file.name = "dummy"
        file.save()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        with open(file.get_temp_pre_path(), "w") as dummy_pcl:
            dummy_pcl.write("This is a dummy file for tests to operate upon.")

        input_file_path = file.get_temp_pre_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch],
                       "input_file_path": input_file_path}

        job_context = array_express._determine_brainarray_package(job_context)
        failure_reason_template = (
            "Unable to read Affy header in input file "
            "/home/user/data_store/temp/TEST2/AFFY_TO_PCL/batch_{0}/dummy "
            "while running AFFY_TO_PCL due to error: Error in (function "
            '(filename, info = c("basic", "full"), verbose = FALSE)  : \n  '
            "Is /home/user/data_store/temp/TEST2/AFFY_TO_PCL/batch_{1}/dummy "
            "really a CEL file? tried reading as text, gzipped text, binary, "
            "gzipped binary, command console and gzipped command console "
            "formats\n"
        )

        self.assertFalse("brainarray_package" in job_context)
        self.assertFalse(job_context["success"])
        self.assertFalse(os.path.isfile(input_file_path))
        self.assertEqual(processor_job.failure_reason,
                         failure_reason_template.format(batch.id, batch.id))


class RunScanUPCTestCase(TestCase):
    def test_success(self):
        batch = init_objects()
        file = batch.files[0]
        batch.platform_accession_code = "TEST"
        batch.save()
        file.internal_location = "TEST/AFFY_TO_PCL"
        file.name = "GSM1426186_UC_colon_inactive_201.CEL"
        file.save()

        # Prevent the test file/directory from getting removed.
        file.remove_temp_directory = MagicMock()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # We have a test file in the repo, but it needs to be in the
        # correct location which depends on the ID of the Batch, which
        # changes based on the order tests are run in.
        test_file_path = "/home/user/data_store/temp/TEST/AFFY_TO_PCL/" + file.name
        input_file_path = file.get_temp_pre_path()
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        shutil.copyfile(test_file_path, input_file_path)

        output_file_path = file.get_temp_post_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch],
                       "brainarray_package": "hugene10sthsentrezgprobe",
                       "input_file_path": input_file_path,
                       "output_file_path": output_file_path}

        # If output_file exists, remove it first.
        if os.path.isfile(output_file_path):
            os.remove(output_file_path)

        job_context = array_express._run_scan_upc(job_context)

        # success is only populated by this function on an error
        self.assertFalse("success" in job_context)
        self.assertTrue(os.path.isfile(output_file_path))

        # Clean up the processed file
        os.remove(output_file_path)

        # Clean up the copied file
        os.remove(input_file_path)

    def test_pcl_concordance(self):
        """Tests that the PCL file generated by SCAN.UPC is valid.

        Runs SCAN.UPC on a sample CEL file, then compares the
        resulting PCL file against a reference PCL file. There are two
        comparisons performed: one which tests that the two files have
        a correlation greater than .99 and another that makes sure
        that the median difference between the reference file and the
        generated file is less than 10% of the standard deviation.
        """
        batch = init_objects()
        # The job operates on the only file associated with the batch.
        file = batch.files[0]
        batch.platform_accession_code = "TEST3"
        batch.save()
        file.internal_location = "TEST3/AFFY_TO_PCL"
        file.name = "GSM955680_DNA10204-001.CEL"
        file.save()

        # Prevent the test file/directory from getting removed.
        file.remove_temp_directory = MagicMock()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        # We have a test file in the repo, but it needs to be in the
        # correct location which depends on the ID of the Batch, which
        # changes based on the order tests are run in.
        test_file_directory = "/home/user/data_store/temp/TEST3/AFFY_TO_PCL/"
        test_file_path = test_file_directory + file.name
        input_file_path = file.get_temp_pre_path()
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        shutil.copyfile(test_file_path, input_file_path)

        output_file_path = file.get_temp_post_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch],
                       "input_file_path": input_file_path,
                       "output_file_path": output_file_path}

        job_context = array_express._determine_brainarray_package(job_context)

        # If output_file exists, remove it first.
        if os.path.isfile(output_file_path):
            os.remove(output_file_path)

        job_context = array_express._run_scan_upc(job_context)

        # success is only populated by this function on an error
        self.assertFalse("success" in job_context)
        self.assertTrue(os.path.isfile(output_file_path))

        # Verification of the PCL file
        reference_file = test_file_directory + "E-GEOD-39088_reference.pcl"

        with closing(open(output_file_path)) as f:
            reference_sample = f.readline().replace("\n", "")

        reference_matrix = pd.read_csv(reference_file, delimiter="\t")
        reference_matrix.rename(index=str,
                                columns={"Unnamed: 0": "genes",
                                         reference_sample: "reference_sample"},
                                inplace=True)

        test_matrix = pd.read_csv(output_file_path, delimiter="\t", skiprows=1, header=None)
        test_matrix.rename(index=str, columns={0: "genes", 1: "test_sample"}, inplace=True)

        intersection = pd.merge(reference_matrix, test_matrix, how="inner", on="genes")
        intersection["epsilon"] = intersection["reference_sample"] - intersection["test_sample"]

        standard_deviations = intersection.std(axis=0)
        min_standard_deviation = min(standard_deviations["reference_sample"],
                                     standard_deviations["test_sample"])
        error_margin = min_standard_deviation * 0.1

        median_epsilon = intersection["epsilon"].abs().median()
        self.assertLess(median_epsilon, error_margin)

        correlation = intersection[["reference_sample", "test_sample"]].corr()
        self.assertGreater(correlation["reference_sample"]["test_sample"], 0.99)

        # Clean up the processed file
        os.remove(output_file_path)

        # Clean up the copied file
        os.remove(input_file_path)

    def test_failure(self):
        batch = init_objects()
        file = batch.files[0]
        batch.platform_accession_code = "TEST2"
        batch.save()
        file.internal_location = "TEST2/AFFY_TO_PCL"
        file.name = "dummy"
        file.save()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        with open(file.get_temp_pre_path(), "w") as dummy_pcl:
            dummy_pcl.write("This is a dummy file for tests to operate upon.")

        input_file_path = file.get_temp_pre_path()
        output_file_path = file.get_temp_post_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch],
                       "brainarray_package": "hugene10sthsentrezgprobe",
                       "input_file_path": input_file_path,
                       "output_file_path": output_file_path}

        job_context = array_express._run_scan_upc(job_context)
        failure_reason_template = (
            'Encountered error in R code while running AFFY_TO_PCL pipeline '
            'during processing of '
            '/home/user/data_store/temp/TEST2/AFFY_TO_PCL/batch_{0}/dummy: '
            'Error in {{ : \n  task 1 failed - '
            '"Is /home/user/data_store/temp/TEST2/AFFY_TO_PCL/batch_{1}/dummy '
            'really a CEL file? tried reading as text, gzipped text, binary, '
            'gzipped binary, command console and gzipped command console '
            'formats\n"\n'
        )

        self.assertFalse(job_context["success"])
        self.assertFalse(os.path.isfile(input_file_path))
        self.assertFalse(os.path.isfile(output_file_path))
        self.assertEqual(processor_job.failure_reason,
                         failure_reason_template.format(batch.id, batch.id))

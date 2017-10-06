import os
from django.test import TestCase
from unittest.mock import MagicMock, patch
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import no_op


def init_objects():
    survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
    survey_job.save()

    batch = Batch(
        survey_job=survey_job,
        source_type="ARRAY_EXPRESS",
        pipeline_required="NO_OP",
        platform_accession_code="TEST",
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
                raw_format="txt",
                processed_format="txt",
                name="test.txt",
                internal_location="TEST/NO_OP/",
                batch=batch)
    file.save()

    batch.files = [file]

    return batch


class RunNoOpTestCase(TestCase):
    @patch("data_refinery_common.models.File.objects")
    def test_success(self, mock_file_objects):
        batch = init_objects()

        # Prevent the test file from getting removed.
        batch.files[0].remove_raw_files = MagicMock()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        output_file_path = batch.files[0].get_processed_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch]}

        # If output_file exists, remove it first.
        if os.path.isfile(output_file_path):
            os.remove(output_file_path)

        job_context = no_op._no_op_processor_fn(job_context)

        # success is only populated by this function on an error
        self.assertTrue(job_context["success"])
        self.assertTrue(os.path.isfile(output_file_path))

        # Clean up the processed file
        os.remove(output_file_path)

    def test_failure(self):
        batch = init_objects()
        batch.platform_accession_code = "TEST2"
        batch.save()
        file = batch.files[0]
        file.internal_location = "TEST2/AFFY_TO_PCL"
        file.name = "dummy"
        file.save()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        output_file_path = file.get_processed_path()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": [batch]}

        # If output_file exists, remove it first.
        if os.path.isfile(output_file_path):
            os.remove(output_file_path)

        job_context = no_op._no_op_processor_fn(job_context)

        # success is only populated by this function on an error
        self.assertFalse(job_context["success"])
        self.assertFalse(os.path.isfile(output_file_path))
        self.assertEqual(processor_job.failure_reason,
                         "Exception caught while moving file dummy")

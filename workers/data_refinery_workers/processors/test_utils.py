import copy
import os
from unittest.mock import MagicMock
from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    File,
    ProcessorJob,
)
from data_refinery_workers.processors import utils


def init_objects():
    survey_job = SurveyJob(
        source_type="ARRAY_EXPRESS"
    )
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
    batch2 = copy.deepcopy(batch)
    batch.save()
    batch2.save()

    file = File(
        size_in_bytes=0,
        download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip/GSM1426072_CD_colon_active_2.CEL",  # noqa
        raw_format="CEL",
        processed_format="PCL",
        name="CE1234.CEL",
        internal_location="A-AFFY-1/AFFY_TO_PCL/",
        batch=batch
    )
    file2 = copy.deepcopy(file)
    file2.name = "CE2345.CEL"
    file2.batch = batch2
    file.save()
    file2.save()

    batch.files = [file]
    batch2.files = [file2]

    return (batch, batch2)


class StartJobTestCase(TestCase):
    def test_success(self):
        batch, batch2 = init_objects()

        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch, batch2])

        job_context = utils.start_job({"job": processor_job})
        # start_job preserves the "job" key
        self.assertEqual(job_context["job"], processor_job)

        # start_job finds the batches and returns them
        self.assertEqual(len(job_context["batches"]), 2)

    def test_failure(self):
        """Fails because there are no batches for the job."""
        processor_job = ProcessorJob()
        processor_job.save()

        job_context = utils.start_job({"job": processor_job})
        self.assertFalse(job_context["success"])


class EndJobTestCase(TestCase):
    def test_success(self):
        batch, batch2 = init_objects()

        processor_job = ProcessorJob()
        processor_job.save()

        utils.end_job({"job": processor_job,
                       "batches": [batch, batch2]})

        processor_job.refresh_from_db()
        self.assertTrue(processor_job.success)
        self.assertIsNotNone(processor_job.end_time)

        batches = Batch.objects.all()
        for batch in batches:
            self.assertEqual(batch.status, BatchStatuses.PROCESSED.value)

    def test_failure(self):
        batch, batch2 = init_objects()

        processor_job = ProcessorJob()
        processor_job.save()

        utils.end_job({"success": False,
                       "job": processor_job,
                       "batches": [batch, batch2]})

        processor_job.refresh_from_db()
        self.assertFalse(processor_job.success)
        self.assertIsNotNone(processor_job.end_time)

        batches = Batch.objects.all()
        for batch in batches:
            self.assertEqual(batch.status, BatchStatuses.DOWNLOADED.value)


class UploadProcessedFilesTestCase(TestCase):
    def setUp(self):
        self.batch, _ = init_objects()

    def tearDown(self):
        expected_path = "/home/user/data_store/processed/A-AFFY-1/AFFY_TO_PCL/CE1234.PCL"
        if os.path.isfile(expected_path):
            os.remove(expected_path)

    def test_success(self):
        file = self.batch.files[0]
        processor_job = ProcessorJob.create_job_and_relationships(batches=[self.batch])
        os.makedirs(file.get_temp_dir(), exist_ok=True)
        with open(file.get_temp_post_path(), "w") as dummy_pcl:
            dummy_pcl.write("This is a dummy file for tests to operate upon.")

        # Verify file was created correctly or else the test which
        # verifies that it was removed won't actually be testing
        # anything
        self.assertTrue(os.path.isfile(file.get_temp_post_path()))

        job_context = {"batches": [self.batch],
                       "job": processor_job,
                       "job_id": processor_job.id}
        job_context = utils.upload_processed_files(job_context)

        self.assertFalse("success" in job_context)
        self.assertTrue(os.path.isfile(file.get_processed_path()))
        self.assertFalse(os.path.isfile(file.get_temp_post_path()))

    def test_failure(self):
        file = self.batch.files[0]
        processor_job = ProcessorJob.create_job_and_relationships(batches=[self.batch])
        job_context = {"batches": [self.batch],
                       "job": processor_job,
                       "job_id": processor_job.id}
        job_context = utils.upload_processed_files(job_context)

        self.assertFalse(job_context["success"])
        self.assertEqual(type(job_context["job"].failure_reason), str)
        self.assertFalse(os.path.isfile(file.get_processed_path()))


class RunPipelineTestCase(TestCase):
    def test_no_job(self):
        mock_processor = MagicMock()
        utils.run_pipeline({"job_id": 100}, [mock_processor])
        mock_processor.assert_not_called()

    def test_processor_failure(self):
        processor_job = ProcessorJob()
        processor_job.save()
        job_context = {"job_id": processor_job.id,
                       "job": processor_job,
                       "batches": []}

        mock_processor = MagicMock()
        mock_processor.__name__ = "Fake processor."
        return_context = copy.copy(job_context)
        return_context["success"] = False
        mock_processor.return_value = return_context

        utils.run_pipeline(job_context, [mock_processor])
        self.assertEqual(mock_processor.call_count, 1)
        processor_job.refresh_from_db()
        self.assertFalse(processor_job.success)
        self.assertIsNotNone(processor_job.end_time)

    def test_value_passing(self):
        """The keys added to job_context and returned by processors will be
        passed through to other processors.
        """
        batch, _ = init_objects()
        processor_job = ProcessorJob.create_job_and_relationships(batches=[batch])

        mock_processor = MagicMock()
        mock_context = {"something_to_pass_along": True,
                        "job": processor_job,
                        "batches": [batch]}
        mock_processor.return_value = mock_context

        def processor_function(job_context):
            self.assertTrue(job_context["something_to_pass_along"])
            return job_context

        test_processor = MagicMock(side_effect=processor_function)

        utils.run_pipeline({"job_id": processor_job.id},
                           [utils.start_job,
                            mock_processor,
                            test_processor,
                            utils.end_job])

        processor_job.refresh_from_db()
        self.assertTrue(processor_job.success)
        self.assertIsNotNone(processor_job.end_time)

        batch.refresh_from_db()
        self.assertEqual(batch.status, BatchStatuses.PROCESSED.value)

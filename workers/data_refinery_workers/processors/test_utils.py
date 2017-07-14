import copy
from unittest.mock import patch, MagicMock
from django.test import TestCase
from data_refinery_models.models import (
    SurveyJob,
    Batch,
    BatchStatuses,
    DownloaderJob,
    DownloaderJobsToBatches,
    ProcessorJob,
    ProcessorJobsToBatches
)
from data_refinery_workers.processors import utils


def init_batch():
    survey_job = SurveyJob(
        source_type="ARRAY_EXPRESS"
    )
    survey_job.save()

    return Batch(
        survey_job=survey_job,
        source_type="ARRAY_EXPRESS",
        size_in_bytes=0,
        download_url="ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip/GSM1426072_CD_colon_active_2.CEL",  # noqa
        raw_format="CEL",
        processed_format="PCL",
        pipeline_required="AFFY_TO_PCL",
        platform_accession_code="A-AFFY-1",
        experiment_accession_code="E-MTAB-3050",
        experiment_title="It doesn't really matter.",
        name="CE1234.CEL",
        internal_location="A-AFFY-1/MICRO_ARRAY_TO_PCL/",
        organism_id=9606,
        organism_name="HOMO SAPIENS",
        release_date="2017-05-05",
        last_uploaded_date="2017-05-05",
        status=BatchStatuses.DOWNLOADED.value
    )


class StartJobTestCase(TestCase):
    def test_success(self):
        batch = init_batch()
        batch2 = copy.deepcopy(batch)
        batch2.name = "CE2345.CEL"
        batch.save()
        batch2.save()

        processor_job = ProcessorJob()
        processor_job.save()
        processor_job_to_batch = ProcessorJobsToBatches(batch=batch,
                                                        processor_job=processor_job)
        processor_job_to_batch.save()
        processor_job_to_batch2 = ProcessorJobsToBatches(batch=batch2,
                                                         processor_job=processor_job)
        processor_job_to_batch2.save()

        kwargs = utils.start_job({"job": processor_job})
        # start_job preserves the "job" key
        self.assertEqual(kwargs["job"], processor_job)

        # start_job finds the batches and returns them
        self.assertEqual(len(kwargs["batches"]), 2)

    def test_failure(self):
        """Fails because there are no batches for the job."""
        processor_job = ProcessorJob()
        processor_job.save()

        kwargs = utils.start_job({"job": processor_job})
        self.assertFalse(kwargs["success"])


class EndJobTestCase(TestCase):
    def test_success(self):
        batch = init_batch()
        batch2 = copy.deepcopy(batch)
        batch2.name = "CE2345.CEL"
        batch.save()
        batch2.save()

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
        batch = init_batch()
        batch2 = copy.deepcopy(batch)
        batch2.name = "CE2345.CEL"
        batch.save()
        batch2.save()

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


class RunPipelineTestCase(TestCase):
    def test_no_job(self):
        mock_processor = MagicMock()
        utils.run_pipeline({"job_id": 100}, [mock_processor])
        mock_processor.assert_not_called()

    def test_processor_failure(self):
        processor_job = ProcessorJob()
        processor_job.save()
        job_dict = {"job_id": processor_job.id,
                    "job": processor_job}

        mock_processor = MagicMock()
        mock_processor.__name__ = "Fake processor."
        return_dict = copy.copy(job_dict)
        return_dict["success"] = False
        mock_processor.return_value = return_dict

        utils.run_pipeline(job_dict, [mock_processor])
        self.assertEqual(mock_processor.call_count, 1)
        processor_job.refresh_from_db()
        self.assertFalse(processor_job.success)
        self.assertIsNotNone(processor_job.end_time)

    def test_value_passing(self):
        """The keys added to kwargs and returned by processors will be
        passed through to other processors.
        """
        batch = init_batch()
        batch.save()
        processor_job = ProcessorJob()
        processor_job.save()
        processor_jobs_to_batches = ProcessorJobsToBatches(batch=batch,
                                                           processor_job=processor_job)
        processor_jobs_to_batches.save()

        mock_processor = MagicMock()
        mock_dict = {"something_to_pass_along": True,
                     "job": processor_job,
                     "batches": [batch]}
        mock_processor.return_value = mock_dict

        def processor_function(kwargs):
            self.assertTrue(kwargs["something_to_pass_along"])
            return kwargs

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

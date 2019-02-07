import copy
import os
from io import StringIO
from unittest.mock import MagicMock
from django.core.management import call_command
from django.test import TestCase
from data_refinery_common.models import (
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    SurveyJob,
    ProcessorJob,
    Sample,
    OriginalFileSampleAssociation,
    Organism
)
from django.utils import timezone
from data_refinery_workers.processors import utils


def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
    pj.save()

    original_file = OriginalFile()
    original_file.source_filename = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"
    original_file.filename = "GSM1426071_CD_colon_active_1.CEL"
    original_file.absolute_file_path = "/home/user/data_store/raw/TEST/CEL/GSM1426071_CD_colon_active_1.CEL"
    original_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = original_file
    assoc1.processor_job = pj
    assoc1.save()

    c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")
    
    sample = Sample()
    sample.title = "Heyo"
    sample.organism = c_elegans
    sample.is_processed = False
    sample.save()

    ogsa = OriginalFileSampleAssociation()
    ogsa.sample = sample
    ogsa.original_file = original_file
    ogsa.save()

    return pj


class StartJobTestCase(TestCase):
    def test_success(self):

        processor_job = prepare_job()

        job_context = utils.start_job({"job": processor_job})

        # start_job preserves the "job" key
        self.assertEqual(job_context["job"], processor_job)

        job_context['success'] = True
        job_context = utils.end_job(job_context)
        for sample in job_context['samples']:
            self.assertTrue(sample.is_processed)

    def test_failure(self):
        """Fails because there are no files for the job."""
        processor_job = ProcessorJob()
        processor_job.save()

        job_context = utils.start_job({"job": processor_job})
        self.assertFalse(job_context["success"])

    def test_bad_restart(self):

        with self.settings(RUNNING_IN_CLOUD=True):
            job = ProcessorJob()
            job.start_time = timezone.now()
            job.success = True
            job.save()
            job_context = utils.start_job({"job": job})

            job = ProcessorJob()
            job.start_time = timezone.now()
            job.success = False
            job.save()
            job_context = utils.start_job({"job": job})

            self.assertRaises(utils.start_job({"job": job}))

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

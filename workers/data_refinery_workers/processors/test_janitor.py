import csv
import json
import os
import shutil
import sys
import zipfile

from io import StringIO
from unittest.mock import patch, call, MagicMock
from nomad import Nomad
from nomad.api.exceptions import URLNotFoundNomadException

from django.core.management import call_command
from django.test import TestCase, tag
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    ComputationalResult,
    ComputedFile,
    Experiment,
    Organism,
    Sample,
    SampleAnnotation,
    SampleResultAssociation,
    ExperimentSampleAssociation,
    Dataset,
    ProcessorJobDatasetAssociation,
    SampleComputedFileAssociation,
    ComputationalResultAnnotation
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import janitor

LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
JOBS = 10

def prepare_job():

    # Create 10 job directories
    for i in range(0, JOBS):
        os.makedirs(LOCAL_ROOT_DIR + '/processor_job_' + str(i), exist_ok=True)

    # Create a job out of the range with index in it to make sure we
    # don't delete index directories since that's where transcriptome
    # indices get downloaded to.
    os.makedirs(LOCAL_ROOT_DIR + '/processor_job_' + str(JOBS+1) + '_index', exist_ok=True)

    # Save two jobs so that we trigger two special circumstances, one
    # where the job is still running and the other where querying
    # nomad raises an exception.
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.nomad_job_id = "running_job"
    pj.save()

    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.nomad_job_id = "missing_job"
    pj.save()

    pj = ProcessorJob()
    pj.pipeline_applied = "JANITOR"
    pj.save()

    return pj

class JanitorTestCase(TestCase):

    @tag("janitor")
    @patch('data_refinery_workers.processors.janitor.Nomad')
    def test_janitor(self, mock_nomad):
        """ Main tester. """
        def mock_get_job(job_id: str):
            if job_id == "running_job":
                return {"Status": "running"}
            else:
                return {"Status": "dead"}

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = mock_get_job
            return ret_value

        mock_nomad.side_effect = mock_init_nomad
        job = prepare_job()
        final_context = janitor.run_janitor(job.pk)

        for i in range(0, JOBS):
            # The job with id 1 should appear running.
            if i == 1:
                self.assertTrue(os.path.exists(LOCAL_ROOT_DIR + '/processor_job_' + str(i)))
            else:
                self.assertFalse(os.path.exists(LOCAL_ROOT_DIR + '/processor_job_' + str(i)))

        self.assertTrue(os.path.exists(LOCAL_ROOT_DIR + '/processor_job_11_index'))

        # Deleted all the working directories except for the one that's still running.
        self.assertEqual(len(final_context['deleted_items']), JOBS-1)

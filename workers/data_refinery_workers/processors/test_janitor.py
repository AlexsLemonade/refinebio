import csv
import json
import os
import shutil
import sys
import zipfile

from io import StringIO
from unittest.mock import patch, call
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

    for i in range(0, JOBS):
        os.makedirs(LOCAL_ROOT_DIR + '/' + 'processor_job_' + str(i), exist_ok=True)

    pj = ProcessorJob()
    pj.pipeline_applied = "JANITOR"
    pj.save()

    return pj

class JanitorTestCase(TestCase):

    @tag("janitor")
    @patch.object(Nomad, 'job')
    def test_janitor(self, mock_get_job):
        """ Main tester. """
        job = prepare_job()
        final_context = janitor.run_janitor(job.pk)
        
        for i in range(0, JOBS):
            self.assertFalse(os.path.exists(LOCAL_ROOT_DIR + '/' + 'processor_job_' + str(i)))

        self.assertEqual(len(final_context['deleted_items']), 10)

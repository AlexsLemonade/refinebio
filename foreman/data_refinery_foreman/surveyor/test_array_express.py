import json
import datetime
from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism,
    Sample
)
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor

class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value="E-MTAB-3050")
        key_value_pair.save()

        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(name="HOMO_SAPIENS",
                            taxonomy_id=9606,
                            is_scientific_name=True)
        organism.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    @patch('data_refinery_foreman.surveyor.external_source.send_job')
    def test_survey(self, mock_send_task):
        """A Simple test of the ArrayExpress surveyor.
        """
        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        ae_surveyor.survey()

        samples = Sample.objects.all()
        downloader_jobs = DownloaderJob.objects.all()

        # We are expecting this to discoever 5 samples.
        self.assertEqual(samples.count(), 5)
        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 1)

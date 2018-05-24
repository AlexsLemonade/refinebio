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
from data_refinery_foreman.surveyor.geo import (
    GeoSurveyor
)


class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="GEP")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="experiment_accession_code",
                                           value="GSE35186")
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
    def test_geo_survey(self, mock_send_task):
        """ Run the GEO surveyor and make sure we get some files to DL!
        """

        ae_surveyor = GeoSurveyor(self.survey_job)
        ae_surveyor.survey()

        self.assertEqual(124, Sample.objects.all().count())
        downloader_jobs = DownloaderJob.objects.all()
        
        # 124 samples + 2 metadata files that get picked up
        self.assertEqual(126, downloader_jobs.count())

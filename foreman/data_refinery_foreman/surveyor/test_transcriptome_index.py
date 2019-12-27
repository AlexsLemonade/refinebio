import os
import json
import urllib
from unittest.mock import Mock, patch, call
from unittest import skip
from django.test import TestCase
from urllib.request import URLError
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism,
)
from data_refinery_foreman.surveyor.transcriptome_index import (
    TranscriptomeIndexSurveyor,
)


class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblPlants"
        )
        key_value_pair.save()

    @patch("data_refinery_foreman.surveyor.external_source.message_queue.send_job")
    def test_survey(self, mock_send_job):
        surveyor = TranscriptomeIndexSurveyor(self.survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertGreater(downloader_jobs.count(), 50)
        send_job_calls = []
        for downloader_job in downloader_jobs:
            send_job_calls.append(call(Downloaders.TRANSCRIPTOME_INDEX, downloader_job))

        mock_send_job.assert_has_calls(send_job_calls)

    def test_correct_index_location(self):
        """ Tests that the files returned actually exist.

        Uses an organism in the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="Ensembl"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Danio rerio"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(self.survey_job)
        files = surveyor.discover_species()[0]

        for file in files:
            urllib.request.urlopen(file.source_url)

    def test_correct_index_location_metazoa(self):
        """ Tests that the files returned actually exist.

        Tests the Metazoa division instead of the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblMetazoa"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Octopus bimaculoides"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(self.survey_job)
        files = surveyor.discover_species()[0]

        for file in files:
            urllib.request.urlopen(file.source_url)

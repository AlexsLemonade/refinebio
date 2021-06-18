import urllib
from unittest.mock import call, patch

from django.test import TestCase

import vcr

from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import DownloaderJob, Organism, SurveyJob, SurveyJobKeyValue
from data_refinery_foreman.surveyor.transcriptome_index import TranscriptomeIndexSurveyor


class SurveyTestCase(TestCase):
    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.transcriptome.survey.yaml")
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey(self, mock_send_job):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblPlants"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertGreater(downloader_jobs.count(), 50)
        send_job_calls = []
        for downloader_job in downloader_jobs:
            send_job_calls.append(call(Downloaders.TRANSCRIPTOME_INDEX, downloader_job))

        mock_send_job.assert_has_calls(send_job_calls)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.transcriptome.correct_index_location.yaml"
    )
    def test_correct_index_location(self):
        """ Tests that the files returned actually exist.

        Uses an organism in the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="Ensembl"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Danio rerio"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        files = surveyor.discover_species()[0]

        # Make sure the organism object got created by making sure
        # this doesn't raise an exception.
        Organism.objects.get(name="DANIO_RERIO")

        for file in files:
            urllib.request.urlopen(file.source_url)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.transcriptome.correct_location_metazoa.yaml"
    )
    def test_correct_index_location_metazoa(self):
        """ Tests that the files returned actually exist.

        Tests the Metazoa division instead of the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblMetazoa"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Octopus bimaculoides"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        files = surveyor.discover_species()[0]

        for file in files:
            urllib.request.urlopen(file.source_url)

        # Make sure the organism object got created by making sure
        # this doesn't raise an exception.
        Organism.objects.get(name="OCTOPUS_BIMACULOIDES")

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.transcriptome.single_plant.yaml")
    def test_single_plant(self):
        """ Tests that the files returned actually exist.

        Tests the Metazoa division instead of the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblPlants"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Arabidopsis thaliana"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        files = surveyor.discover_species()[0]

        for file in files:
            urllib.request.urlopen(file.source_url)

        # Make sure the organism object got created by making sure
        # this doesn't raise an exception.
        Organism.objects.get(name="ARABIDOPSIS_THALIANA")

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.transcriptome.correct_location_protist.yaml"
    )
    def test_correct_index_location_protist(self):
        """ Tests that the files returned actually exist.

        Tests the Metazoa division instead of the main division.
        """
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblProtists"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="Leishmania major"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        files = surveyor.discover_species()[0]

        for file in files:
            urllib.request.urlopen(file.source_url)

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.transcriptome.survey_fungi.yaml")
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey_fungi(self, mock_send_job):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblFungi"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="CANDIDA_ALBICANS"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertEqual(downloader_jobs.count(), 1)
        send_job_calls = []
        for downloader_job in downloader_jobs:
            send_job_calls.append(call(Downloaders.TRANSCRIPTOME_INDEX, downloader_job))

        mock_send_job.assert_has_calls(send_job_calls)

        # Make sure the organism object got created with the correct
        # taxonomy id by making sure this doesn't raise an exception.
        Organism.objects.get(name="CANDIDA_ALBICANS", taxonomy_id=5476)

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.transcriptome.survey_bacteria.yaml")
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey_bacteria(self, mock_send_job):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblBacteria"
        )
        key_value_pair.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="organism_name", value="PSEUDOMONAS_AERUGINOSA"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertEqual(downloader_jobs.count(), 1)
        send_job_calls = []
        for downloader_job in downloader_jobs:
            send_job_calls.append(call(Downloaders.TRANSCRIPTOME_INDEX, downloader_job))

        mock_send_job.assert_has_calls(send_job_calls)

        # Make sure the organism object got created with the correct
        # taxonomy id by making sure this doesn't raise an exception.
        Organism.objects.get(name="PSEUDOMONAS_AERUGINOSA", taxonomy_id=287)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.transcriptome.survey_bacteria_none.yaml"
    )
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey_bacteria_none(self, mock_send_job):
        """When surveying fungi an organism_name must be supplied."""
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblBacteria"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertEqual(downloader_jobs.count(), 0)

        mock_send_job.assert_not_called()

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.transcriptome.survey_fungi_none.yaml"
    )
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey_fungi_none(self, mock_send_job):
        """When surveying fungi an organism_name must be supplied."""
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="ensembl_division", value="EnsemblFungi"
        )
        key_value_pair.save()

        surveyor = TranscriptomeIndexSurveyor(survey_job)
        surveyor.survey(source_type="TRANSCRIPTOME_INDEX")

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertEqual(downloader_jobs.count(), 0)

        mock_send_job.assert_not_called()

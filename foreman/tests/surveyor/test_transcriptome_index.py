from unittest.mock import call, patch

from django.test import TestCase

import vcr

from data_refinery_common.enums import Downloaders
from data_refinery_common.models import DownloaderJob, Organism, SurveyJob, SurveyJobKeyValue
from data_refinery_foreman.surveyor.transcriptome_index import TranscriptomeIndexSurveyor
from data_refinery_foreman.surveyor.utils import requests_has_content_length


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

    def test_all_correct_index_location(self):
        """
        Tests one organism from all supported Ensembl divisions.
        Fungi and Bacteria organisms have mappings defined
        in config/organism_strain_mapping.csv
        """
        test_organisms = [
            {"division": "Ensembl", "name": "Danio rerio"},
            {"division": "EnsemblMetazoa", "name": "Octopus bimaculoides"},
            {"division": "EnsemblPlants", "name": "Arabidopsis thaliana"},
            {"division": "EnsemblProtists", "name": "Leishmania major"},
            {"division": "EnsemblFungi", "name": "Candida albicans"},
            {"division": "EnsemblBacteria", "name": "Pseudomonas aeruginosa"},
        ]

        for organism in test_organisms:
            scientific_name = organism["name"].upper().replace(" ", "_")
            # Break up cassettes into multiple files.
            with vcr.use_cassette(
                f"/home/user/data_store/cassettes/surveyor.transcriptome.correct_index_location_{scientific_name.lower()}.yaml"
            ):
                survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
                survey_job.save()
                key_value_pair = SurveyJobKeyValue(
                    survey_job=survey_job, key="ensembl_division", value=organism["division"]
                )
                key_value_pair.save()

                key_value_pair = SurveyJobKeyValue(
                    survey_job=survey_job, key="organism_name", value=organism["name"]
                )
                key_value_pair.save()

                surveyor = TranscriptomeIndexSurveyor(survey_job)
                files = surveyor.discover_species()[0]
                # Make sure the organism object got created by making sure
                # this doesn't raise an exception.
                Organism.objects.get(name=scientific_name)

                # Assert that both Fasta and GTF files are found.
                self.assertEqual(len(files), 2)

                # Ensure that the file exists at the location by performing
                # a HEAD request and reading the content-length header.
                checks = [requests_has_content_length(file.source_url) for file in files]
                self.assertTrue(all(checks))

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

        # Make sure the organism object was created with the correct
        # species taxonomy id by making sure this doesn't raise an exception.
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

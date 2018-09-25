import os
import json
from unittest.mock import Mock, patch, call
from django.test import TestCase
from urllib.request import URLError
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import (
    Batch,
    DownloaderJob,
    SurveyJob,
    SurveyJobKeyValue,
    Organism
)
from data_refinery_foreman.surveyor.transcriptome_index import TranscriptomeIndexSurveyor


class SurveyTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="ensembl_division",
                                           value="EnsemblPlants")
        key_value_pair.save()

    @patch('data_refinery_foreman.surveyor.external_source.send_job')
    @patch("data_refinery_foreman.surveyor.transcriptome_index.urllib.request.urlopen")
    @patch("data_refinery_foreman.surveyor.transcriptome_index.requests.get")
    def test_survey(self, mock_get, mock_urlopen, mock_send_job):
        json_file_path = os.path.join(os.path.dirname(__file__), "test_transcriptome_species.json")
        with open(json_file_path, "r") as json_file:
            species_json = json.load(json_file)

        # Insert the organisms into the database so the model doesn't call the
        # taxonomy API to populate them.
        for species in species_json:
            # Account for the subtle difference between the API for
            # the main Ensembl division and the API for the rest of
            # them.
            name_key = "common_name" if "common_name" in species else "name"
            taxonomy_key = "taxonomy_id" if "taxonomy_id" in species else "taxon_id"
            organism = Organism(name=species[name_key].upper(),
                                taxonomy_id=species[taxonomy_key],
                                is_scientific_name=True)
            organism.save()

        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.json.return_value = species_json

        # There are two possible file locations. The correct one is
        # determined by making a request to one to see if it
        # exists. This URLError simulates it not existing.
        mock_urlopen.side_effect = URLError("404 or something")

        surveyor = TranscriptomeIndexSurveyor(self.survey_job)
        surveyor.survey()

        downloader_jobs = DownloaderJob.objects.order_by("id").all()
        self.assertEqual(downloader_jobs.count(), len(species_json))
        send_job_calls = []
        for downloader_job in downloader_jobs:
            send_job_calls.append(
                call(Downloaders.TRANSCRIPTOME_INDEX,
                     downloader_job.id))

        mock_send_job.assert_has_calls(send_job_calls)

        # There should be 2 Batches for each species (long and short
        # transcriptome lengths).
        batches = Batch.objects.all()
        self.assertEqual(batches.count(), len(species_json) * 2)
        # And each batch has two files: fasta and gtf
        for batch in batches:
            self.assertEqual(len(batch.files), 2)

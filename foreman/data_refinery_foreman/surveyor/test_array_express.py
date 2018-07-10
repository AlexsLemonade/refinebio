import json
import datetime
import requests

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

    @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_survey(self, mock_send_task):
        """A Simple test of the ArrayExpress surveyor.
        """
        ae_surveyor = ArrayExpressSurveyor(self.survey_job)
        ae_surveyor.survey()

        samples = Sample.objects.all()
        downloader_jobs = DownloaderJob.objects.all()

        # We are expecting this to discover 5 samples.
        self.assertEqual(samples.count(), 5)
        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 1)

        sample = Sample.objects.first()
        self.assertTrue(' (hgu95av2)' in sample.pretty_platform())

    def test_determine_accession(self):
        """Test of the `determine_sample_accession` function
        """
        ae_surveyor = ArrayExpressSurveyor(self.survey_job)

        EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
        SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"

        ex_accessions = [
                            "E-MTAB-3050",
                            "E-MEXP-669",
                            "E-MEXP-2215",
                            "E-MEXP-2288",
                            "E-MEXP-2381",
                            "E-MTAB-6739",
                        ]

        for ex_accession in ex_accessions:
            samples_endpoint = SAMPLES_URL.format(ex_accession)
            r = requests.get(samples_endpoint, timeout=60)
            samples = r.json()["experiment"]["sample"]

            # An experiment can have many samples
            for sample in samples:

                # For some reason, this sample has no files associated with it.
                if "file" not in sample or len(sample['file']) == 0:
                    continue

                # The accession code is not a simple matter to determine.
                sample_source_name = sample["source"].get("name", "")
                sample_assay_name = sample["assay"].get("name", "")

                has_raw = False
                for sub_file in sample['file']:

                    # For ex: E-GEOD-15645
                    if isinstance(sub_file['comment'], list):
                        sub_file_mod = sub_file
                        sub_file_mod['comment'] = sub_file['comment'][0]
                    else:
                        sub_file_mod = sub_file

                    if sub_file_mod['type'] == "data" and sub_file_mod['comment'].get('value', None) != None:
                        has_raw = True
                    if 'raw' in sub_file_mod['comment'].get('value', ''):
                        has_raw = True

                    # Skip derived data if we have it raw.
                    if has_raw and "derived data" in sub_file['type']:
                        continue
                    elif (not has_raw) and "derived data" not in sub_file['type']:
                        # If there is a platform warning then we don't want raw data.
                        has_raw = False
                        continue
                    filename = sub_file["name"]

                sample_accession_code = ae_surveyor.determine_sample_accession(
                    ex_accession,
                    sample_source_name,
                    sample_assay_name,
                    filename)
                self.assertTrue(sample_accession_code is not None)

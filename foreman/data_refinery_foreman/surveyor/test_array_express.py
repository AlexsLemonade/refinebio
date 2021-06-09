import datetime
from unittest.mock import patch

from django.test import TestCase
from django.utils import timezone

import requests
import vcr

from data_refinery_common.models import (
    DownloaderJob,
    Experiment,
    Organism,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor


class SurveyTestCase(TestCase):
    def setUp(self):
        # Insert human organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    def create_job_for_accession(self, accession_code: str):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value=accession_code
        )
        key_value_pair.save()

        return survey_job

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.array_express.survey.yaml")
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey(self, mock_send_task):
        """A Simple test of the ArrayExpress surveyor."""
        survey_job = self.create_job_for_accession("E-MTAB-3050")
        ae_surveyor = ArrayExpressSurveyor(survey_job)
        ae_surveyor.survey()

        samples = Sample.objects.all()
        downloader_jobs = DownloaderJob.objects.all()

        # We are expecting this to discover 5 samples.
        self.assertEqual(samples.count(), 5)

        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 1)

        experiment = Experiment.objects.first()
        self.assertEqual(experiment.accession_code, "E-MTAB-3050")
        self.assertEqual(
            experiment.source_first_published, datetime.datetime(2014, 10, 31, tzinfo=timezone.utc)
        )
        self.assertEqual(
            experiment.source_last_modified, datetime.datetime(2014, 10, 30, tzinfo=timezone.utc)
        )

        sample = Sample.objects.first()
        self.assertTrue(" (hgu95av2)" in sample.pretty_platform)
        # Confirm the sample's protocol_info
        self.assertEqual(len(sample.protocol_info), 9)
        self.assertEqual(sample.protocol_info[0]["Accession"], "P-MTAB-41854")
        self.assertEqual(sample.protocol_info[0]["Text"], "Aliquoting of biomaterials.")
        self.assertEqual(sample.protocol_info[0]["Type"], "split")

        survey_job2 = self.create_job_for_accession("E-GEOD-44719")
        ae_surveyor = ArrayExpressSurveyor(survey_job2)
        ae_surveyor.survey()

        # We are expecting this to discover 77 samples.
        self.assertEqual(samples.count(), 77 + 5)

        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 2)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.array_express.survey_with_protocol_list.yaml"
    )
    def test_survey_with_protocol_list(self):
        """Tests an edge case that came up after months:
        https://github.com/AlexsLemonade/refinebio/issues/761
        """
        survey_job = self.create_job_for_accession("E-MEXP-2381")
        ae_surveyor = ArrayExpressSurveyor(survey_job)
        ae_surveyor.survey()

        samples = Sample.objects.all()
        downloader_jobs = DownloaderJob.objects.all()

        # We are expecting this to discover 2 samples.
        self.assertEqual(samples.count(), 2)

        # And for one DownloaderJob to be created for all of them.
        self.assertEqual(downloader_jobs.count(), 1)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.array_express.determine_accession.yaml"
    )
    def test_determine_accession(self):
        """Test of the `determine_sample_accession` function
        """
        survey_job = self.create_job_for_accession("E-MTAB-3050")
        ae_surveyor = ArrayExpressSurveyor(survey_job)

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
                if "file" not in sample or len(sample["file"]) == 0:
                    continue

                # The accession code is not a simple matter to determine.
                sample_source_name = sample["source"].get("name", "")
                sample_assay_name = sample["assay"].get("name", "")

                has_raw = False
                for sub_file in sample["file"]:

                    # For ex: E-GEOD-15645
                    if isinstance(sub_file["comment"], list):
                        sub_file_mod = sub_file
                        sub_file_mod["comment"] = sub_file["comment"][0]
                    else:
                        sub_file_mod = sub_file

                    if (
                        sub_file_mod["type"] == "data"
                        and sub_file_mod["comment"].get("value", None) is not None
                    ):
                        has_raw = True
                    if "raw" in sub_file_mod["comment"].get("value", ""):
                        has_raw = True

                    # Skip derived data if we have it raw.
                    if has_raw and "derived data" in sub_file["type"]:
                        continue
                    elif (not has_raw) and "derived data" not in sub_file["type"]:
                        # If there is a platform warning then we don't want raw data.
                        has_raw = False
                        continue
                    filename = sub_file["name"]

                sample_accession_code = ae_surveyor.determine_sample_accession(
                    ex_accession, sample_source_name, sample_assay_name, filename
                )
                self.assertTrue(sample_accession_code is not None)

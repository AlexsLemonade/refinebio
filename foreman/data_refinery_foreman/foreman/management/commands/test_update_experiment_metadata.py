from django.test import TransactionTestCase
from unittest.mock import patch

from data_refinery_common.models import (
    Experiment,
    ExperimentSampleAssociation,
    Sample,
)
from data_refinery_foreman.foreman.management.commands.update_experiment_metadata import Command
from data_refinery_foreman.surveyor.test_sra import mocked_requests_get


class SurveyTestCase(TransactionTestCase):
    def tearDown(self):
        Experiment.objects.all().delete()

    @patch("data_refinery_foreman.surveyor.sra.requests.get")
    def test_sra_experiment_missing_metadata(self, mock_get):
        """Tests that an SRA experiment has its missing metadata added."""
        mock_get.side_effect = mocked_requests_get

        # 1. Create an experiment with a bad title
        BAD_TITLE = "GEO accession GSE1337 is currently private\
 and is scheduled to be released on Jan 01, 1970."

        experiment = Experiment()
        experiment.accession_code = "DRP003977"
        experiment.source_database = "SRA"
        experiment.title = BAD_TITLE
        experiment.save()

        # 2. We need to add a sample because the way that the SRA surveyor finds metadata is
        # through run accessions
        sample = Sample()
        sample.accession_code = "DRR002116"
        sample.technology = "RNA-SEQ"
        sample.source_database = "SRA"
        sample.title = "Not important"
        sample.save()

        ExperimentSampleAssociation.objects.get_or_create(experiment=experiment, sample=sample)

        # 3. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test that the title was fixed
        self.assertNotEqual(
            Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0].title,
            BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes
        command = Command()
        command.handle()

    def test_geo_experiment_missing_metadata(self):
        """Tests that a GEO experiment has its missing metadata added."""

        # 1. Create an experiment with a bad title
        BAD_TITLE = "GEO accession GSE1337 is currently private\
 and is scheduled to be released on Jan 01, 1970."

        experiment = Experiment()
        experiment.accession_code = "GSE11915"
        experiment.source_database = "GEO"
        experiment.title = BAD_TITLE
        experiment.save()

        # 2. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test that the title was fixed
        self.assertNotEqual(
            Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0].title,
            BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes
        command = Command()
        command.handle()

    def test_array_express_experiment_missing_metadata(self):
        """Tests that an ArrayExpress experiment has its missing metadata added."""

        # 1. Create an experiment with a bad title
        BAD_TITLE = "GEO accession GSE1337 is currently private\
 and is scheduled to be released on Jan 01, 1970."

        experiment = Experiment()
        experiment.accession_code = "E-MTAB-3050"
        experiment.source_database = "ARRAY_EXPRESS"
        experiment.title = BAD_TITLE
        experiment.save()

        # 2. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test that the title was fixed
        self.assertNotEqual(
            Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0].title,
            BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes
        command = Command()
        command.handle()

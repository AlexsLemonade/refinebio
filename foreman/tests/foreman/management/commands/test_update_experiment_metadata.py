from django.test import TransactionTestCase

from data_refinery_common.models import Experiment, ExperimentSampleAssociation, Sample
from data_refinery_foreman.foreman.management.commands.update_experiment_metadata import Command


class SurveyTestCase(TransactionTestCase):
    BAD_TITLE = (
        "BAD TITLE: Accession is currently private "
        "and is scheduled to be released on Jan 01, 1970."
    )

    def tearDown(self):
        Experiment.objects.all().delete()

    def test_sra_experiment_missing_metadata(self):
        """Tests that an SRA experiment has its missing metadata added."""

        # 1. Create an experiment with a bad title.
        experiment = Experiment()
        experiment.accession_code = "DRP003977"
        experiment.source_database = "SRA"
        experiment.title = SurveyTestCase.BAD_TITLE
        experiment.save()

        # 2. We need to add a sample because the way that the SRA surveyor finds metadata is
        # through run accessions.
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

        # Test that the title was fixed.
        experiment = Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0]
        self.assertNotEqual(
            experiment.title,
            SurveyTestCase.BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes.
        command = Command()
        command.handle()

    def test_sra_experiment_missing_alternate_accession(self):
        """Tests that an SRA experiment has its missing alternate_accession_code added."""

        # 1. Create an experiment without an alternate_accession_code.
        experiment = Experiment()
        experiment.accession_code = "SRP094947"
        experiment.source_database = "SRA"
        experiment.title = "Not important"
        experiment.save()

        # 2. We need to add a sample because the way that the SRA surveyor finds metadata is
        # through run accessions.
        sample = Sample()
        sample.accession_code = "SRR5099111"
        sample.technology = "RNA-SEQ"
        sample.source_database = "SRA"
        sample.title = "Not important"
        sample.save()

        ExperimentSampleAssociation.objects.get_or_create(experiment=experiment, sample=sample)

        # 3. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # 4. Refresh the experiment
        experiment.refresh_from_db()

        # Test that the correct alternate_accession_code was added
        self.assertEquals(experiment.alternate_accession_code, "GSE92260")

    def test_geo_experiment_missing_metadata(self):
        """Tests that a GEO experiment has its missing metadata added."""

        # 1. Create an experiment with a bad title.
        experiment = Experiment()
        experiment.accession_code = "GSE11915"
        experiment.source_database = "GEO"
        experiment.title = SurveyTestCase.BAD_TITLE
        experiment.save()

        # 2. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test that the title was fixed.
        experiment = Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0]
        self.assertNotEqual(
            experiment.title,
            SurveyTestCase.BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes.
        command = Command()
        command.handle()

    def test_array_express_experiment_missing_metadata(self):
        """Tests that an ArrayExpress experiment has its missing metadata added."""

        # 1. Create an experiment with a bad title
        experiment = Experiment()
        experiment.accession_code = "E-MTAB-3050"
        experiment.source_database = "ARRAY_EXPRESS"
        experiment.title = SurveyTestCase.BAD_TITLE
        experiment.save()

        # 2. Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test that the title was fixed.
        experiment = Experiment.objects.get_or_create(accession_code=experiment.accession_code)[0]
        self.assertNotEqual(
            experiment.title,
            SurveyTestCase.BAD_TITLE,
        )

        # Run the command again to make sure that it does not fail if there are no changes.
        command = Command()
        command.handle()

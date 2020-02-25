from django.core.management import call_command
from django.test import TestCase

from data_refinery_common.models import (
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.foreman.management.commands.rerun_salmon_old_samples import (
    update_salmon_all_experiments,
)


def setup_experiments() -> None:
    """Creates three experiments for testing purposes.

    One experiment will not have a GSE* accession code.
    Both experiments with GSE* accession codes will have two
    samples. One has a sample that has incorrect platform information
    so the experiment will need to be re-surveyed.
    """
    organism = Organism.get_object_for_name("HOMO_SAPIENS")

    # Experiment that needs to be re-surveyed
    experiment = Experiment.objects.create(
        accession_code="GSE12417", technology="MICROARRAY", source_database="GEO"
    )

    # Correct platform
    sample = Sample.objects.create(
        accession_code="GSM311750",
        source_database="GEO",
        technology="MICROARRAY",
        platform_accession_code="hgu133a",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    # Incorrect Platform
    sample = Sample.objects.create(
        accession_code="GSM316652",
        organism=organism,
        source_database="GEO",
        technology="MICROARRAY",
        platform_accession_code="hgu133a",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    # Experiment that does not need to be re-surveyed
    experiment = Experiment.objects.create(
        accession_code="GSE9890", technology="MICROARRAY", source_database="GEO"
    )

    sample = Sample.objects.create(
        accession_code="GSM249671",
        source_database="GEO",
        technology="MICROARRAY",
        platform_accession_code="hgu133plus2",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    sample = Sample.objects.create(
        accession_code="GSM249672",
        organism=organism,
        source_database="GEO",
        technology="MICROARRAY",
        platform_accession_code="hgu133plus2",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    # Experiment that isn't even the right source_database.
    experiment = Experiment.objects.create(
        accession_code="SRP12345", technology="RNA-SEQ", source_database="SRA"
    )

    sample = Sample.objects.create(
        accession_code="SRR123145",
        organism=organism,
        source_database="SRA",
        technology="RNA-SEQ",
        platform_accession_code="IlluminaHiSeq1000",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    # This second sample is just to make checking things easy because
    # all experiments start with 2 samples in these tests.
    sample = Sample.objects.create(
        accession_code="SRR123146",
        organism=organism,
        source_database="SRA",
        technology="RNA-SEQ",
        platform_accession_code="IlluminaHiSeq1000",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)


class CorrectAffyCdfTestCase(TestCase):
    """
    Tests that new processor jobs are created for samples that belong to experiments that were
    processed with multiple versions of Salmon
    """

    def test_dry_run(self):
        setup_experiments()
        # Setup is done, actually run the command.
        args = []
        options = {"dry_run": True}
        call_command("correct_affy_cdfs", *args, **options)

        # The experiments actually have more than 2 samples each, so
        # if any of them have more than 2 samples then it got
        # resurveyed and the dry run did not work.
        for experiment in Experiment.objects.all():
            self.assertEqual(2, experiment.samples.count())

    def test_resurveying(self):
        setup_experiments()
        # Setup is done, actually run the command.
        args = []
        options = {}
        call_command("correct_affy_cdfs", *args, **options)

        # Verify that the two experiments with correct platform information do not get queued.
        unchanged_experiments = Experiment.objects.filter(
            accession_code__in=["GSE9890", "SRP12345"]
        )

        for experiment in unchanged_experiments:
            self.assertEqual(2, experiment.samples.count())

        # GSE12417 should have been unsurveyed.
        with self.assertRaises(Experiment.DoesNotExist):
            Experiment.objects.get(accession_code="GSE12417")

        survey_jobs = SurveyJob.objects.all()
        self.assertEqual(1, survey_jobs.count())

        self.assertEqual("GSE12417", SurveyJobKeyValue.objects.first().value)

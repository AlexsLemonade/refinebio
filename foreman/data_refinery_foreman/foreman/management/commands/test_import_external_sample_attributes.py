from unittest.mock import patch

from django.test import TestCase

import vcr

from data_refinery_common.models import (
    Contribution,
    Experiment,
    ExperimentSampleAssociation,
    OntologyTerm,
    Sample,
    SampleAttribute,
)
from data_refinery_foreman.foreman.management.commands.import_external_sample_attributes import (
    Command,
    import_metadata,
    import_sample_attributes,
)

TEST_METADATA = "/home/user/data_store/externally_supplied_metadata/test_data/metadata.json"


class ImportExternalSampleAttributesTestCase(TestCase):
    def setUp(self):
        experiment = Experiment()
        experiment.accession_code = "GSE000"
        experiment.alternate_accession_code = "E-GEOD-000"
        experiment.title = "NONONONO"
        experiment.description = "Boooooourns. Wasabi."
        experiment.technology = "RNA-SEQ"
        experiment.save()
        self.experiment = experiment

        # Create some samples to attach metadata to
        sample = Sample()
        sample.accession_code = "SRR123"
        sample.technology = "RNA-SEQ"
        sample.source_database = "SRA"
        sample.title = "Not important"
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        sample2 = Sample()
        sample2.accession_code = "SRR456"
        sample2.technology = "RNA-SEQ"
        sample2.source_database = "SRA"
        sample2.title = "Not important"
        sample2.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample2
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        # Create the ontology terms I'm using in the tests
        name = OntologyTerm()
        name.ontology_term = "PATO:0000122"
        name.human_readable_name = "length"
        name.save()

        unit = OntologyTerm()
        unit.ontology_term = "UO:0010012"
        unit.human_readable_name = "thou"
        unit.save()

        contribution = Contribution()
        contribution.source_name = "refinebio_tests"
        contribution.methods_url = "ccdatalab.org"
        contribution.save()
        self.contribution = contribution

    #
    # Test import_sample_attributes()
    #

    def test_skip_unknown_sample(self):
        """Make sure that if someone has metadata for a sample that we haven't
        surveyed then we just do nothing"""

        METADATA = [{"PATO:0000122": {"value": 25, "unit": "UO:0010012"}}]
        import_sample_attributes("SRR789", METADATA, self.contribution)
        self.assertEqual(SampleAttribute.objects.all().count(), 0)

    def test_import_invalid_ontology_term(self):
        METADATA = [{"PATO:0000122": {"value": 25, "unit": "thou"}}]
        self.assertRaises(
            ValueError, import_sample_attributes, "SRR123", METADATA, self.contribution
        )

        METADATA = [{"length": {"value": 25, "unit": "UO:0010012"}}]
        self.assertRaises(
            ValueError, import_sample_attributes, "SRR123", METADATA, self.contribution
        )

    def test_import_valid_sample_attributes(self):
        METADATA = [{"PATO:0000122": {"value": 25, "unit": "UO:0010012"}}]
        import_sample_attributes("SRR123", METADATA, self.contribution)

        self.assertEqual(SampleAttribute.objects.all().count(), 1)

        contributed_metadata = Sample.objects.get(accession_code="SRR123").contributed_metadata
        self.assertEqual(
            contributed_metadata[self.contribution.source_name]["length"],
            {"unit": "thou", "value": 25},
        )

    #
    # Test import_metadata()
    #

    def test_import_valid_metadata(self):
        METADATA = [
            {
                "sample_accession": "SRR123",
                "attributes": [{"PATO:0000122": {"value": 25, "unit": "UO:0010012"}}],
            }
        ]

        import_metadata(METADATA, self.contribution)

        self.assertEqual(SampleAttribute.objects.all().count(), 1)

        contributed_metadata = Sample.objects.get(accession_code="SRR123").contributed_metadata
        self.assertEqual(
            contributed_metadata[self.contribution.source_name]["length"],
            {"unit": "thou", "value": 25},
        )

    #
    # End-to-end test
    #
    @vcr.use_cassette("/home/user/data_store/cassettes/foreman.sample_attributes.end-to-end.yaml")
    def test_management_command(self):
        sample = Sample()
        sample.accession_code = "DRR001173"
        sample.technology = "RNA-SEQ"
        sample.source_database = "SRA"
        sample.title = "Not important"
        sample.save()

        command = Command()
        SOURCE_NAME = "refinebio_tests"
        command.handle(file=TEST_METADATA, source_name=SOURCE_NAME, methods_url="ccdatalab.org")

        self.assertEqual(SampleAttribute.objects.all().count(), 1)

        contributed_metadata = sample.contributed_metadata
        self.assertEqual(
            set(contributed_metadata[SOURCE_NAME]["biological sex"].keys()),
            {"value", "confidence"},
        )
        self.assertEqual(
            contributed_metadata[SOURCE_NAME]["biological sex"]["value"].human_readable_name,
            "female",
        )

        self.assertAlmostEqual(
            contributed_metadata[SOURCE_NAME]["biological sex"]["confidence"], 0.7856624891880539
        )

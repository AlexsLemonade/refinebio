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

        metadata = Sample.objects.get(accession_code="SRR123").to_metadata_dict()
        self.assertIsNotNone(metadata.get("other_metadata", None))
        self.assertEqual(len(metadata["other_metadata"]), 1)
        self.assertEqual(metadata["other_metadata"][0]["name"]["term"], "PATO:0000122")
        self.assertEqual(metadata["other_metadata"][0]["name"]["name"], "length")
        self.assertEqual(metadata["other_metadata"][0]["unit"]["term"], "UO:0010012")
        self.assertEqual(metadata["other_metadata"][0]["unit"]["name"], "thou")
        self.assertEqual(metadata["other_metadata"][0]["value"], 25)

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

        metadata = Sample.objects.get(accession_code="SRR123").to_metadata_dict()
        self.assertIsNotNone(metadata.get("other_metadata", None))
        self.assertEqual(len(metadata["other_metadata"]), 1)
        self.assertEqual(metadata["other_metadata"][0]["name"]["term"], "PATO:0000122")
        self.assertEqual(metadata["other_metadata"][0]["name"]["name"], "length")
        self.assertEqual(metadata["other_metadata"][0]["unit"]["term"], "UO:0010012")
        self.assertEqual(metadata["other_metadata"][0]["unit"]["name"], "thou")
        self.assertEqual(metadata["other_metadata"][0]["value"], 25)
        self.assertEqual(metadata["other_metadata"][0]["probability"], "unknown")

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
        command.handle(
            file=TEST_METADATA, source_name="refinebio_tests", methods_url="ccdatalab.org"
        )

        self.assertEqual(SampleAttribute.objects.all().count(), 1)

        metadata = sample.to_metadata_dict()
        self.assertIsNotNone(metadata.get("other_metadata", None))
        self.assertEqual(len(metadata["other_metadata"]), 1)
        # Make sure everything matches what was in TEST_METADATA
        self.assertEqual(metadata["other_metadata"][0]["name"]["term"], "PATO:0000047")
        self.assertEqual(metadata["other_metadata"][0]["name"]["name"], "biological sex")
        self.assertEqual(metadata["other_metadata"][0]["value"]["term"], "PATO:0000383")
        self.assertEqual(metadata["other_metadata"][0]["value"]["name"], "female")
        self.assertAlmostEqual(metadata["other_metadata"][0]["probability"], 0.7856624891880539)

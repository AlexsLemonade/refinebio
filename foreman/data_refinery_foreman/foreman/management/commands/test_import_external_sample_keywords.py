from unittest.mock import patch

from django.test import TestCase

import vcr

from data_refinery_common.models import (
    Experiment,
    ExperimentSampleAssociation,
    OntologyTerm,
    Sample,
    SampleKeyword,
)
from data_refinery_foreman.foreman.management.commands.import_external_sample_keywords import (
    Command,
    import_keywords,
)

SUBMITTER = "Refinebio tests"
TEST_KEYWORDS = "/home/user/data_store/externally_supplied_metadata/test_data/keywords.json"


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

        # Create some samples to attach keywords to
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

    def test_skip_unknown_sample(self):
        """Make sure that if someone has keywords for a sample that we haven't
        surveyed then we just do nothing"""

        KEYWORDS = {"SRR789": ["PATO:0000122", "UO:0010012"]}
        import_keywords(KEYWORDS, SUBMITTER)
        self.assertEqual(SampleKeyword.objects.all().count(), 0)

    def test_import_invalid_ontology_term(self):
        KEYWORDS = {"SRR123": ["width", "UO:0010012"]}
        self.assertRaises(ValueError, import_keywords, KEYWORDS, SUBMITTER)

    def test_import_valid_sample_keywords(self):
        KEYWORDS = {"SRR123": ["PATO:0000122", "UO:0010012"]}
        import_keywords(KEYWORDS, SUBMITTER)

        self.assertEqual(SampleKeyword.objects.all().count(), 2)

        sample = Sample.objects.get(accession_code="SRR123")
        self.assertEqual(
            set(sample.keywords.values_list("name__human_readable_name", flat=True)),
            set(["length", "thou"]),
        )

    #
    # End-to-end test
    #
    @vcr.use_cassette("/home/user/data_store/cassettes/foreman.sample_keywords.end-to-end.yaml")
    def test_management_command(self):
        sample = Sample()
        sample.accession_code = "DRR000897"
        sample.technology = "RNA-SEQ"
        sample.source_database = "SRA"
        sample.title = "Not important"
        sample.save()

        sample2 = Sample()
        sample2.accession_code = "DRR001173"
        sample2.technology = "RNA-SEQ"
        sample2.source_database = "SRA"
        sample2.title = "Not important"
        sample2.save()

        command = Command()
        command.handle(file=TEST_KEYWORDS, submitter=SUBMITTER)

        # If you look below you'll only see 14, but this is because DRR001173
        # has two pairs of terms from different ontologies with the same
        # human-readable name
        self.assertEqual(SampleKeyword.objects.all().count(), 16)

        # I checked all of these manually in the Ontology Lookup Service, the
        # data itself comes from MetaSRA

        self.assertEqual(
            set(sample.keywords.values_list("name__human_readable_name", flat=True)),
            set(["late embryonic stage", "serum", "late embryo", "cultured cell"]),
        )

        self.assertEqual(
            set(sample2.keywords.values_list("name__human_readable_name", flat=True)),
            set(
                [
                    "epithelial neoplasm",
                    "cancer",
                    "neoplasm",
                    "bladder disease",
                    "bladder carcinoma",
                    "carcinoma",
                    "disease",
                    "disease of cellular proliferation",
                    "cell line",
                    "cultured cell",
                ]
            ),
        )

from unittest.mock import patch

from django.test import TestCase

import vcr

from data_refinery_common.models.ontology_term import OntologyTerm, get_human_readable_name_from_api


class TestOntologyTerm(TestCase):
    @vcr.use_cassette("/home/user/data_store/cassettes/common.ontology_term.import_from_api.yaml")
    def test_get_or_create_from_api_fetches(self):
        self.assertEqual(
            "month", OntologyTerm.get_or_create_from_api("UO:0000035").human_readable_name
        )

    @vcr.use_cassette("/home/user/data_store/cassettes/common.ontology_term.invalid.yaml")
    def test_get_or_create_invalid_term(self):
        self.assertRaises(ValueError, OntologyTerm.get_or_create_from_api, "UO:9999999")

    @patch("data_refinery_common.models.ontology_term.get_human_readable_name_from_api")
    def test_get_or_create_from_api_returns_cached(self, mock_api_call):
        mock_api_call.return_value = "medulloblastoma"

        self.assertEqual(
            "medulloblastoma",
            OntologyTerm.get_or_create_from_api("EFO:0002939").human_readable_name,
        )

        mock_api_call.assert_called_once_with("EFO:0002939")
        mock_api_call.reset_mock()

        self.assertEqual(
            "medulloblastoma",
            OntologyTerm.get_or_create_from_api("EFO:0002939").human_readable_name,
        )

        mock_api_call.assert_not_called()

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/common.ontology_term.import_entire_ontology.yaml"
    )
    @patch("data_refinery_common.models.ontology_term.get_human_readable_name_from_api")
    def test_import_entire_ontology(self, mock_api_call):
        # We shouldn't be hitting the API at all here, because we should have
        # the ontology already imported
        mock_api_call.return_value = "The wrong answer"

        # I used the UO ontology here because it is much smaller than other important
        # ontologies like EFO, which could take upwards of a minute to download and parse
        created_terms = OntologyTerm.import_entire_ontology("uo")

        self.assertEqual(
            OntologyTerm.objects.all().count(),
            559,  # Since we are using a VCR, this number should not change until we refresh it
        )
        self.assertEqual(OntologyTerm.objects.all().count(), created_terms)

        self.assertEqual(
            "month", OntologyTerm.get_or_create_from_api("UO:0000035").human_readable_name
        )

        mock_api_call.assert_not_called()

    @vcr.use_cassette("/home/user/data_store/cassettes/common.ontology_term.import_cl.yaml")
    @patch("data_refinery_common.models.ontology_term.get_human_readable_name_from_api")
    def test_import_cl(self, mock_api_call):
        """Try importing the CL ontology, which was updated at some point and
        broke the original parsing code"""

        # We shouldn't be hitting the API at all here, because we should have
        # the ontology already imported
        mock_api_call.return_value = "The wrong answer"

        created_terms = OntologyTerm.import_entire_ontology("cl")

        self.assertEqual(
            OntologyTerm.objects.all().count(),
            2493,  # Since we are using a VCR, this number should not change until we refresh it
        )
        self.assertEqual(OntologyTerm.objects.all().count(), created_terms)

        self.assertEqual(
            "hematopoietic cell",
            OntologyTerm.get_or_create_from_api("CL:0000988").human_readable_name,
        )

        mock_api_call.assert_not_called()

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/common.ontology_term.import_cellosaurus.yaml"
    )
    @patch("data_refinery_common.models.ontology_term.get_human_readable_name_from_api")
    def test_import_cellosaurus(self, mock_api_call):
        """The cellosaurus ontology is not part of the OLS so we need to handle
        it separately. We still include it because metaSRA uses its terms.


        NOTE: the actual cellosaurus ontology is massive, so this VCR was
        created using a trimmed-down cellosaurus ontology where I went in and
        deleted a bunch of the publications and cell lines from the respective
        lists. Besides this, no alterations were made to the original."""

        # We shouldn't be hitting the API at all here, because we should have
        # the ontology already imported
        mock_api_call.return_value = "The wrong answer"

        created_terms = OntologyTerm.import_entire_ontology("cvcl")

        self.assertEqual(
            OntologyTerm.objects.all().count(),
            34,  # This is the number I counted in the file
        )
        self.assertEqual(OntologyTerm.objects.all().count(), created_terms)

        self.assertEqual(
            "#W7079",
            OntologyTerm.get_or_create_from_api("CVCL:E549").human_readable_name,
        )

        mock_api_call.assert_not_called()

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/common.ontology_term.cellosaurus_import_from_api.yaml"
    )
    def test_get_or_create_from_cellosaurus_api(self):
        self.assertEqual(
            "LNCaP", OntologyTerm.get_or_create_from_api("CVCL:0395").human_readable_name
        )

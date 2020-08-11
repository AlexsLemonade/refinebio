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

        term = OntologyTerm()
        term.ontology_term = "EFO:0002939"
        term.human_readable_name = "medulloblastoma"
        term.save()

        self.assertEqual(
            "medulloblastoma",
            OntologyTerm.get_or_create_from_api("EFO:0002939").human_readable_name,
        )

        mock_api_call.assert_not_called()

    @patch("data_refinery_common.models.ontology_term.get_human_readable_name_from_api")
    def test_poke_term(self, mock_api_call):
        mock_api_call.return_value = "medulloblastoma"

        OntologyTerm.poke_term("EFO:0002939")

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
        OntologyTerm.import_entire_ontology("uo")

        self.assertEqual(
            "month", OntologyTerm.get_or_create_from_api("UO:0000035").human_readable_name
        )

        mock_api_call.assert_not_called()

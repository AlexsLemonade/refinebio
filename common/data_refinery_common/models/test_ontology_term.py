from unittest.mock import patch
from django.test import TestCase
from data_refinery_common.models.ontology_term import (
    OntologyTerm, get_human_readable_name_from_api,
)


def mocked_get_human_readable_name(ontology_term):
    if ontology_term == "EFO:0002939":
        return "medulloblastoma"
    elif ontology_term == "UO:0000035":
        return "month"
    else:
        raise Exception("Was not expecting the term: " + ontology_term)


class TestOntologyTerm(TestCase):
    def tearDown(self):
        OntologyTerm.objects.all().delete()

    def test_get_human_readable_name(self):
        self.assertEqual("medulloblastoma",
                         get_human_readable_name_from_api("EFO:0002939"))
        self.assertEqual("month",
                         get_human_readable_name_from_api("UO:0000035"))

    @patch('data_refinery_common.models.ontology_term.get_human_readable_name_from_api')
    def test_get_or_create_fetches(self, mock_api_call):
        mock_api_call.side_effect = mocked_get_human_readable_name

        self.assertEqual("month", OntologyTerm.get_or_create("UO:0000035").human_readable_name)

        mock_api_call.assert_called_once_with("UO:0000035")

    @patch('data_refinery_common.models.ontology_term.get_human_readable_name_from_api')
    def test_get_or_create_returns_cached(self, mock_api_call):
        mock_api_call.side_effect = mocked_get_human_readable_name

        term = OntologyTerm()
        term.ontology_term = "EFO:0002939"
        term.human_readable_name = "medulloblastoma"
        term.save()

        self.assertEqual("medulloblastoma",
                         OntologyTerm.get_or_create("EFO:0002939").human_readable_name)

        mock_api_call.assert_not_called()

    @patch('data_refinery_common.models.ontology_term.get_human_readable_name_from_api')
    def test_poke_term(self, mock_api_call):
        mock_api_call.side_effect = mocked_get_human_readable_name

        OntologyTerm.poke_term("EFO:0002939")

        self.assertEqual("medulloblastoma",
                         OntologyTerm.get_or_create("EFO:0002939").human_readable_name)

        mock_api_call.assert_called_once_with("EFO:0002939")

    @patch('data_refinery_common.models.ontology_term.get_human_readable_name_from_api')
    def test_import_entire_ontology(self, mock_api_call):
        # I used the UO ontology here because it is much smaller than other important
        # ontologies like EFO, which could take upwards of a minute to download and parse
        OntologyTerm.import_entire_ontology("uo")

        self.assertEqual("month", OntologyTerm.get_or_create("UO:0000035").human_readable_name)

        mock_api_call.assert_not_called()

from unittest.mock import Mock, call, patch

from django.test import TestCase

from data_refinery_common.models.organism import (
    EFETCH_URL,
    ESEARCH_URL,
    InvalidNCBITaxonomyId,
    Organism,
    UnknownOrganismId,
)

ESEARCH_RESPONSE_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN"
    "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult>
  <Count>1</Count>
  <RetMax>1</RetMax>
  <RetStart>0</RetStart>
  <IdList>
    <Id>9606</Id>
  </IdList>
  <TranslationSet/>
  <TranslationStack>
    <TermSet>
      <Term>homo sapiens[Scientific Name]</Term>
      <Field>Scientific Name</Field>
      <Count>1</Count>
      <Explode>N</Explode>
    </TermSet>
    <OP>GROUP</OP>
  </TranslationStack>
  <QueryTranslation>homo sapiens[Scientific Name]</QueryTranslation>
</eSearchResult>"""

ESEARCH_NOT_FOUND_XML = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN"
    "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">
<eSearchResult>
  <Count>0</Count>
  <RetMax>0</RetMax>
  <RetStart>0</RetStart>
  <IdList/>
  <TranslationSet/>
  <QueryTranslation>(man[Scientific Name])</QueryTranslation>
  <ErrorList>
    <PhraseNotFound>blah</PhraseNotFound>
  </ErrorList>
  <WarningList>
    <OutputMessage>No items found.</OutputMessage>
  </WarningList>
</eSearchResult>"""

EFETCH_RESPONSE_XML = """<?xml version="1.0" ?>
<!DOCTYPE TaxaSet PUBLIC "-//NLM//DTD Taxon, 14th January 2002//EN"
"https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd">
<TaxaSet><Taxon>
    <TaxId>9606</TaxId>
    <ScientificName>Homo sapiens</ScientificName>
    <OtherNames>
        <GenbankCommonName>human</GenbankCommonName>
        <CommonName>man</CommonName>
        <Name>
            <ClassCDE>authority</ClassCDE>
            <DispName>Homo sapiens Linnaeus, 1758</DispName>
        </Name>
    </OtherNames>
    <ParentTaxId>9605</ParentTaxId>
    <Rank>species</Rank>
    <Division>Primates</Division>
    <GeneticCode>
        <GCId>1</GCId>
        <GCName>Standard</GCName>
    </GeneticCode>
    <MitoGeneticCode>
        <MGCId>2</MGCId>
        <MGCName>Vertebrate Mitochondrial</MGCName>
    </MitoGeneticCode>
    <Lineage>cellular organisms</Lineage>
    <LineageEx>
        <Taxon>
            <TaxId>131567</TaxId>
            <ScientificName>cellular organisms</ScientificName>
            <Rank>no rank</Rank>
        </Taxon>
    </LineageEx>
    <CreateDate>1995/02/27 09:24:00</CreateDate>
    <UpdateDate>2017/02/28 16:38:58</UpdateDate>
    <PubDate>1992/05/26 01:00:00</PubDate>
</Taxon>

</TaxaSet>"""

EFETCH_NOT_FOUND_XML = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eEfetchResult PUBLIC "-//NLM//DTD efetch 20131226//EN"
"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20131226/efetch.dtd">
<eFetchResult>
<ERROR>ID list is empty! Possibly it has no correct IDs.</ERROR>
</eFetchResult>"""


def mocked_requests_get(url, parameters):
    mock = Mock(ok=True)
    if url is not ESEARCH_URL:
        mock.text = "This is wrong."
    else:
        try:
            if parameters["field"] == "scin":
                mock.text = ESEARCH_NOT_FOUND_XML
            else:
                mock.text = "This is also wrong."
        except KeyError:
            mock.text = ESEARCH_RESPONSE_XML

    return mock


class OrganismModelTestCase(TestCase):
    def tearDown(self):
        Organism.objects.all().delete()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_cached_names_are_found(self, mock_get):
        Organism.objects.create(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)

        name = Organism.get_name_for_id(9606)

        self.assertEqual(name, "HOMO_SAPIENS")
        mock_get.assert_not_called()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_cached_ids_are_found(self, mock_get):
        Organism.objects.create(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)

        id = Organism.get_id_for_name("Homo Sapiens")

        self.assertEqual(id, 9606)
        mock_get.assert_not_called()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_uncached_scientific_names_are_found(self, mock_get):
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = ESEARCH_RESPONSE_XML

        taxonomy_id = Organism.get_id_for_name("Homo Sapiens")

        self.assertEqual(taxonomy_id, 9606)
        mock_get.assert_called_once_with(
            ESEARCH_URL,
            {
                "db": "taxonomy",
                "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09",
                "field": "scin",
                "term": "HOMO_SAPIENS",
            },
        )

        # The first call should have stored the organism record in the
        # database so this call should not make a request.
        mock_get.reset_mock()
        new_id = Organism.get_id_for_name("Homo Sapiens")

        self.assertEqual(new_id, 9606)
        mock_get.assert_not_called()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_uncached_other_names_are_found(self, mock_get):
        mock_get.side_effect = mocked_requests_get

        taxonomy_id = Organism.get_id_for_name("Human")

        self.assertEqual(taxonomy_id, 9606)
        mock_get.assert_has_calls(
            [
                call(
                    ESEARCH_URL,
                    {
                        "db": "taxonomy",
                        "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09",
                        "field": "scin",
                        "term": "HUMAN",
                    },
                ),
                call(
                    ESEARCH_URL,
                    {
                        "db": "taxonomy",
                        "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09",
                        "term": "HUMAN",
                    },
                ),
            ]
        )

        # The first call should have stored the organism record in the
        # database so this call should not make a request.
        mock_get.reset_mock()
        new_id = Organism.get_id_for_name("Human")

        self.assertEqual(new_id, 9606)
        mock_get.assert_not_called()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_uncached_ids_are_found(self, mock_get):
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = EFETCH_RESPONSE_XML

        organism_name = Organism.get_name_for_id(9606)

        self.assertEqual(organism_name, "HOMO_SAPIENS")
        mock_get.assert_called_once_with(
            EFETCH_URL,
            {"db": "taxonomy", "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09", "id": "9606"},
        )

        # The first call should have stored the organism record in the
        # database so this call should not make a request.
        mock_get.reset_mock()
        new_name = Organism.get_name_for_id(9606)

        self.assertEqual(new_name, "HOMO_SAPIENS")
        mock_get.assert_not_called()

    @patch("data_refinery_common.models.organism.requests.get")
    def test_invalid_ids_cause_exceptions(self, mock_get):
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = EFETCH_NOT_FOUND_XML

        with self.assertRaises(InvalidNCBITaxonomyId):
            Organism.get_name_for_id(0)

    @patch("data_refinery_common.models.organism.requests.get")
    def test_unfound_names_raise(self, mock_get):
        """If we can't find the taxonomy id, it's likely a bad organism name."""
        mock_get.return_value = Mock(ok=True)
        mock_get.return_value.text = ESEARCH_NOT_FOUND_XML

        with self.assertRaises(UnknownOrganismId):
            Organism.get_id_for_name("blah")

        mock_get.assert_has_calls(
            [
                call(
                    ESEARCH_URL,
                    {
                        "db": "taxonomy",
                        "field": "scin",
                        "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09",
                        "term": "BLAH",
                    },
                ),
                call(
                    ESEARCH_URL,
                    {
                        "db": "taxonomy",
                        "api_key": "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09",
                        "term": "BLAH",
                    },
                ),
            ]
        )

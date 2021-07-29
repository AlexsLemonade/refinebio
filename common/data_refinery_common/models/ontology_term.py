import xml.etree.ElementTree as ET

from django.db import models

import requests

from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

# The Ontology Lookup Service URL used to look up a single ontology term
OLS_TERM_URL_TEMPLATE = "http://www.ebi.ac.uk/ols/api/terms?id={}"

# The Ontology Lookup Service
OLS_ONTOLOGY_INFO_URL_TEMPLATE = "https://www.ebi.ac.uk/ols/api/ontologies/{}"

# The Ontology Lookup Service URL used to download an entire ontology
OLS_ONTOLOGY_URL_TEMPLATE = "https://www.ebi.ac.uk/ols/ontologies/{}/download"

CELLOSAURUS_TERM_URL_TEMPLATE = "https://web.expasy.org/cellosaurus/{}.txt"
CELLOSAURUS_ONTOLOGY_URL = "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml"


def get_human_readable_name_from_cellosaurus(ontology_term: str) -> str:
    response = requests.get(CELLOSAURUS_TERM_URL_TEMPLATE.format(ontology_term.replace(":", "_")))

    for line in response.iter_lines():
        line = line.decode("utf-8").lstrip("<pre>").rstrip("</pre>")
        cols = line.split("   ")
        if cols[0] == "ID":
            return cols[1]

    raise ValueError("Could not find a human-readable name for '{}'".format(ontology_term))


def get_human_readable_name_from_ols(ontology_term: str) -> str:
    response = requests.get(OLS_TERM_URL_TEMPLATE.format(ontology_term))

    if response.json().get("_embedded", None) is None:
        raise ValueError("Could not find a human-readable name for '{}'".format(ontology_term))

    terms = response.json()["_embedded"]["terms"]

    # The same term can appear in multiple databases, but we only want the human-readable name so
    # it doesn't matter which database we get it from
    for term in terms:
        # Sanity check if for some reason the first item returned isn't actually the one we want
        if term["obo_id"] == ontology_term:
            return term["label"]

    raise ValueError("We can't find {} in the Ontology Lookup Service".format(ontology_term))


def get_human_readable_name_from_api(ontology_term: str) -> str:
    if get_ontology_prefix(ontology_term) == "cvcl":
        return get_human_readable_name_from_cellosaurus(ontology_term)
    elif is_ols_ontology(get_ontology_prefix(ontology_term)):
        return get_human_readable_name_from_ols(ontology_term)
    else:
        raise ValueError(
            f"Tried to get human readable name for term {ontology_term} in unknown ontology."
        )


def get_ontology_prefix(ontology_term: str) -> str:
    """Returns the ontology prefix for an ontology term

    Example: For EFO:0002939, it would return efo
    """
    return ontology_term.split(":")[0].lower()


def is_ols_ontology(ontology_prefix: str) -> bool:
    return requests.get(OLS_ONTOLOGY_INFO_URL_TEMPLATE.format(ontology_prefix)).status_code == 200


class OntologyTerm(models.Model):
    """The mapping between a human-readable name and an ontology term"""

    class Meta:
        db_table = "ontology_terms"

    ontology_term = models.TextField(unique=True)
    human_readable_name = models.TextField()

    def to_dict(self):
        return {
            "term": self.ontology_term,
            "name": self.human_readable_name,
        }

    @staticmethod
    def import_ols_ontology(ontology_prefix: str) -> int:
        response = requests.get(OLS_ONTOLOGY_URL_TEMPLATE.format(ontology_prefix))
        ontology_xml = ET.fromstring(response.text)

        namespace = {
            "owl": "http://www.w3.org/2002/07/owl#",
            "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
        }

        new_terms = 0

        # For some reason, there are two ways to find ontology terms in OWL files.
        # The first is <owl:Class> tags with an <rdfs:label> child
        for child in ontology_xml.findall("owl:Class/[rdfs:label]", namespace):
            about = child.attrib.get("{" + namespace["rdf"] + "}about")
            ontology_term = about.split("/")[-1].replace("_", ":")
            human_readable_name = child.find("rdfs:label", namespace).text

            if get_ontology_prefix(ontology_term) == ontology_prefix:
                term, created = OntologyTerm.objects.get_or_create(ontology_term=ontology_term)
                term.human_readable_name = human_readable_name
                term.save()

                new_terms += 1 if created else 0

        # The other way is <rdf:Description> tags with an rdf:about attribute
        # and a <rdfs:label> child
        for child in ontology_xml.findall("rdf:Description[@rdf:about]/[rdfs:label]", namespace):
            about = child.attrib.get("{" + namespace["rdf"] + "}about")
            ontology_term = about.split("/")[-1].replace("_", ":")
            human_readable_name = child.find("rdfs:label", namespace).text

            if get_ontology_prefix(ontology_term) == ontology_prefix:
                term, created = OntologyTerm.objects.get_or_create(ontology_term=ontology_term)
                term.human_readable_name = human_readable_name
                term.save()

                new_terms += 1 if created else 0

        return new_terms

    @staticmethod
    def import_cellosaurus_ontology() -> int:
        response = requests.get(CELLOSAURUS_ONTOLOGY_URL)
        ontology_xml = ET.fromstring(response.text)

        new_terms = 0

        for cell_line in ontology_xml.findall(
            # Find all <cell-line>'s that are a child of a <cell-line-list>
            "cell-line-list/cell-line/"
            # that have a primary accession
            "accession-list/accession[@type='primary']/../../"
            # and an identifier name
            "name-list/name[@type='identifier']/../.."
            # (see https://docs.python.org/3/library/xml.etree.elementtree.html#xpath-support)
        ):

            primary_accession = cell_line.find("accession-list/accession[@type='primary']")
            ontology_term = primary_accession.text.replace("_", ":")
            identifier = cell_line.find("name-list/name[@type='identifier']")

            term, created = OntologyTerm.objects.get_or_create(ontology_term=ontology_term)
            term.human_readable_name = identifier.text
            term.save()

            new_terms += 1 if created else 0

        return new_terms

    @classmethod
    def import_entire_ontology(cls, ontology_prefix: str) -> int:
        """Given an ontology prefix, download the entire ontology and import it"""

        # Normalize ontology prefixes to be stored in lowercase
        ontology_prefix = ontology_prefix.lower()

        if ontology_prefix == "cvcl":
            return cls.import_cellosaurus_ontology()
        elif is_ols_ontology(ontology_prefix):
            return cls.import_ols_ontology(ontology_prefix)
        else:
            raise ValueError(
                f"Cannot import the {ontology_prefix} ontology; no source database recognized."
            )

    @staticmethod
    def _create_from_api(ontology_term: str) -> "OntologyTerm":
        """Creates and saves an OntologyTerm instance by querying the OLS API"""
        term = OntologyTerm()
        term.ontology_term = ontology_term
        term.human_readable_name = get_human_readable_name_from_api(ontology_term)
        term.save()
        return term

    @classmethod
    def get_or_create_from_api(cls, ontology_term: str) -> "OntologyTerm":
        match = cls.objects.filter(ontology_term=ontology_term).first()
        if match is not None:
            return match

        return cls._create_from_api(ontology_term)

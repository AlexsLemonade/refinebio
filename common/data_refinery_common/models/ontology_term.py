import xml.etree.ElementTree as ET

from django.db import models

import requests

from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

# The Ontology Lookup Service URL used to look up a single ontology term
TERM_URL_TEMPLATE = "http://www.ebi.ac.uk/ols/api/terms?id={}"

# The Ontology Lookup Service URL used to download an entire ontology
ONTOLOGY_URL_TEMPLATE = "https://www.ebi.ac.uk/ols/ontologies/{}/download"


def get_human_readable_name_from_api(ontology_term: str) -> str:
    response = requests.get(TERM_URL_TEMPLATE.format(ontology_term))

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


class OntologyTerm(models.Model):
    """ The mapping between a human-readable name and an ontology term """

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
    def _get_ontology_prefix(ontology_term: str) -> str:
        """ Returns the ontology prefix for an ontology term

        Example: For EFO:0002939, it would return efo
        """
        return ontology_term.split(":")[0].lower()

    @staticmethod
    def import_entire_ontology(ontology_prefix: str):
        """ Given an ontology prefix, download the entire ontology and import it """
        response = requests.get(ONTOLOGY_URL_TEMPLATE.format(ontology_prefix.lower()))
        ontology_xml = ET.fromstring(response.text)

        namespace = {
            "owl": "http://www.w3.org/2002/07/owl#",
            "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
        }

        # For some reason, there are two ways to find ontology terms in OWL files.
        # The first is <owl:Class> tags with an <rdfs:label> child
        for child in ontology_xml.findall("owl:Class/[rdfs:label]", namespace):
            about = child.attrib.get("{" + namespace["rdf"] + "}about")
            ontology_term = about.split("/")[-1].replace("_", ":")
            human_readable_name = child.find("rdfs:label", namespace).text

            if OntologyTerm._get_ontology_prefix(ontology_term) == ontology_prefix:
                term, _ = OntologyTerm.objects.get_or_create(ontology_term=ontology_term)
                term.human_readable_name = human_readable_name
                term.save()

        # The other way is <rdf:Description> tags with an <rdfs:label> child
        for child in ontology_xml.findall("rdf:Description/[rdfs:label]", namespace):
            about = child.attrib.get("{" + namespace["rdf"] + "}about")

            # Something seems to have changed and this no longer to
            # match anything. I'm leaving it on the off chance that it
            # ever does pick anything up.
            if about:
                ontology_term = about.split("/")[-1].replace("_", ":")
                human_readable_name = child.find("rdfs:label", namespace).text

                if OntologyTerm._get_ontology_prefix(ontology_term) == ontology_prefix:
                    term, _ = OntologyTerm.objects.get_or_create(ontology_term=ontology_term)
                    term.human_readable_name = human_readable_name
                    term.save()

    @staticmethod
    def _create_from_api(ontology_term: str) -> "OntologyTerm":
        """ Creates and saves an OntologyTerm instance by querying the OLS API """
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

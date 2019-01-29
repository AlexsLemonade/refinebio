import requests
from xml.etree import ElementTree

from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimeTrackedModel


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


NCBI_ROOT_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ESEARCH_URL = NCBI_ROOT_URL + "esearch.fcgi"
EFETCH_URL = NCBI_ROOT_URL + "efetch.fcgi"
TAXONOMY_DATABASE = "taxonomy"


class UnscientificNameError(Exception):
    pass


class InvalidNCBITaxonomyId(Exception):
    pass


def get_scientific_name(taxonomy_id: int) -> str:
    parameters = {"db": TAXONOMY_DATABASE, "id": str(taxonomy_id)}
    response = requests.get(EFETCH_URL, parameters)

    root = ElementTree.fromstring(response.text)
    taxon_list = root.findall("Taxon")

    if len(taxon_list) == 0:
        logger.error("No names returned by ncbi.nlm.nih.gov for organism "
                     + "with taxonomy ID %d.",
                     taxonomy_id)
        raise InvalidNCBITaxonomyId

    return taxon_list[0].find("ScientificName").text.upper()


def get_taxonomy_id(organism_name: str) -> int:
    parameters = {"db": TAXONOMY_DATABASE, "term": organism_name}
    response = requests.get(ESEARCH_URL, parameters)

    root = ElementTree.fromstring(response.text)
    id_list = root.find("IdList").findall("Id")

    if len(id_list) == 0:
        logger.error("Unable to retrieve NCBI taxonomy ID number for organism "
                     + "with name: %s",
                     organism_name)
        return 0
    elif len(id_list) > 1:
        logger.warn("Organism with name %s returned multiple NCBI taxonomy ID "
                    + "numbers.",
                    organism_name)

    return int(id_list[0].text)


def get_taxonomy_id_scientific(organism_name: str) -> int:
    parameters = {"db": TAXONOMY_DATABASE, "field": "scin", "term": organism_name}
    response = requests.get(ESEARCH_URL, parameters)

    print(response.text)
    root = ElementTree.fromstring(response.text)
    id_list = root.find("IdList").findall("Id")

    if len(id_list) == 0:
        raise UnscientificNameError
    elif len(id_list) > 1:
        logger.warn("Organism with name %s returned multiple NCBI taxonomy ID "
                    + "numbers.",
                    organism_name)

    return int(id_list[0].text)


class Organism(models.Model):
    """Provides a lookup between organism name and taxonomy ids.

    Should only be used via the two class methods get_name_for_id and
    get_id_for_name. These methods will populate the database table
    with any missing values by accessing the NCBI API.
    """

    def __str__ (self):
        return str(self.name)

    name = models.CharField(max_length=256, unique=True)
    taxonomy_id = models.IntegerField(unique=True)
    is_scientific_name = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(Organism, self).save(*args, **kwargs)

    @classmethod
    def get_name_for_id(cls, taxonomy_id: int) -> str:
        try:
            organism = (cls.objects
                        .filter(taxonomy_id=taxonomy_id)
                        .order_by("-is_scientific_name")
                        [0])
        except IndexError:
            name = get_scientific_name(taxonomy_id).upper().replace(' ', '_')
            organism = Organism(name=name,
                                taxonomy_id=taxonomy_id,
                                is_scientific_name=True)
            organism.save()

        return organism.name

    @classmethod
    def get_id_for_name(cls, name: str) -> id:
        name = name.upper().replace(' ', '_')
        try:
            organism = (cls.objects
                        .filter(name=name)
                        [0])
        except IndexError:
            is_scientific_name = False
            try:
                taxonomy_id = get_taxonomy_id_scientific(name)
                is_scientific_name = True
            except UnscientificNameError:
                taxonomy_id = get_taxonomy_id(name)

            organism = Organism(name=name,
                                taxonomy_id=taxonomy_id,
                                is_scientific_name=is_scientific_name)
            organism.save()

        return organism.taxonomy_id

    @classmethod
    def get_object_for_name(cls, name: str) -> id:
        name = name.upper()
        name = name.replace(' ', '_')
        try:
            organism = (cls.objects
                        .filter(name=name)
                        [0])
        except IndexError:
            is_scientific_name = False
            try:
                taxonomy_id = get_taxonomy_id_scientific(name)
                is_scientific_name = True
            except UnscientificNameError:
                taxonomy_id = get_taxonomy_id(name)

            organism = Organism(name=name,
                                taxonomy_id=taxonomy_id,
                                is_scientific_name=is_scientific_name)
            organism.save()

        return organism


    class Meta:
        db_table = "organisms"

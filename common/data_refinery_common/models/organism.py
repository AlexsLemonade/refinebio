import requests
from xml.etree import ElementTree

from django.apps import apps
from django.db import models
from django.utils import timezone


from data_refinery_common.models.base_models import TimeTrackedModel


from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)


NCBI_ROOT_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = get_env_variable(
    "NCBI_API_KEY", "3a1f8d818b0aa05d1aa3c334fa2cc9a17e09"
)  # This is only used by eUtils and for organisms that aren't cached yet - it's harmless to share.
ESEARCH_URL = NCBI_ROOT_URL + "esearch.fcgi"
EFETCH_URL = NCBI_ROOT_URL + "efetch.fcgi"
TAXONOMY_DATABASE = "taxonomy"


class UnscientificNameError(Exception):
    pass


class InvalidNCBITaxonomyId(Exception):
    pass


def get_scientific_name(taxonomy_id: int) -> str:
    parameters = {
        "db": TAXONOMY_DATABASE,
        "id": str(taxonomy_id),
        "api_key": NCBI_API_KEY,
    }
    response = requests.get(EFETCH_URL, parameters)

    try:
        root = ElementTree.fromstring(response.text)
        taxon_list = root.findall("Taxon")
    except Exception as e:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(taxon_list) == 0:
        logger.error(
            "No names returned by ncbi.nlm.nih.gov for organism "
            + "with taxonomy ID %d.",
            taxonomy_id,
        )
        raise InvalidNCBITaxonomyId

    return taxon_list[0].find("ScientificName").text.upper()


def get_taxonomy_id(organism_name: str) -> int:
    parameters = {
        "db": TAXONOMY_DATABASE,
        "term": organism_name,
        "api_key": NCBI_API_KEY,
    }
    response = requests.get(ESEARCH_URL, parameters)

    try:
        root = ElementTree.fromstring(response.text)
        id_list = root.find("IdList").findall("Id")
    except Exception as e:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(id_list) == 0:
        logger.error(
            "Unable to retrieve NCBI taxonomy ID number for organism "
            + "with name: %s",
            organism_name,
        )
        return 0
    elif len(id_list) > 1:
        logger.warn(
            "Organism with name %s returned multiple NCBI taxonomy ID " + "numbers.",
            organism_name,
        )

    return int(id_list[0].text)


def get_taxonomy_id_scientific(organism_name: str) -> int:
    parameters = {
        "db": TAXONOMY_DATABASE,
        "field": "scin",
        "term": organism_name,
        "api_key": NCBI_API_KEY,
    }
    response = requests.get(ESEARCH_URL, parameters)

    try:
        root = ElementTree.fromstring(response.text)
        id_list = root.find("IdList").findall("Id")
    except Exception as e:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(id_list) == 0:
        raise UnscientificNameError
    elif len(id_list) > 1:
        logger.warn(
            "Organism with name %s returned multiple NCBI taxonomy ID " + "numbers.",
            organism_name,
        )

    return int(id_list[0].text)


class Organism(models.Model):
    """Provides a lookup between organism name and taxonomy ids.

    Should only be used via the two class methods get_name_for_id and
    get_id_for_name. These methods will populate the database table
    with any missing values by accessing the NCBI API.
    """

    class Meta:
        db_table = "organisms"

    name = models.CharField(max_length=256, unique=True)
    taxonomy_id = models.IntegerField(unique=True)
    is_scientific_name = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    experiments = models.ManyToManyField(
        "Experiment", through="ExperimentOrganismAssociation"
    )
    qn_target = models.ForeignKey(
        "ComputationalResult", blank=True, null=True, on_delete=models.SET_NULL
    )

    def __str__(self):
        return str(self.name)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(Organism, self).save(*args, **kwargs)

    def get_genus(self):
        return self.name.split("_")[0]

    @classmethod
    def get_name_for_id(cls, taxonomy_id: int) -> str:
        try:
            organism = cls.objects.filter(taxonomy_id=taxonomy_id).order_by(
                "-is_scientific_name"
            )[0]
        except IndexError:
            name = get_scientific_name(taxonomy_id).upper().replace(" ", "_")
            organism = Organism(
                name=name, taxonomy_id=taxonomy_id, is_scientific_name=True
            )
            organism.save()

        return organism.name

    @classmethod
    def get_id_for_name(cls, name: str) -> id:
        name = name.upper().replace(" ", "_")
        try:
            organism = cls.objects.filter(name=name)[0]
        except IndexError:
            is_scientific_name = False
            try:
                taxonomy_id = get_taxonomy_id_scientific(name)
                is_scientific_name = True
            except UnscientificNameError:
                taxonomy_id = get_taxonomy_id(name)

            organism = Organism(
                name=name,
                taxonomy_id=taxonomy_id,
                is_scientific_name=is_scientific_name,
            )
            organism.save()

        return organism.taxonomy_id

    @classmethod
    def get_object_for_name(cls, name: str):
        name = name.upper()
        name = name.replace(" ", "_")
        try:
            organism = cls.objects.filter(name=name)[0]
        except IndexError:
            is_scientific_name = False
            try:
                taxonomy_id = get_taxonomy_id_scientific(name)
                is_scientific_name = True
            except UnscientificNameError:
                taxonomy_id = get_taxonomy_id(name)

            organism = Organism(
                name=name,
                taxonomy_id=taxonomy_id,
                is_scientific_name=is_scientific_name,
            )
            organism.save()

        return organism

    @classmethod
    def get_objects_with_qn_targets(cls):
        """ Return a list of Organisms who already have valid QN targets associated with them. """
        return Organism.objects.all().filter(qn_target__isnull=False)

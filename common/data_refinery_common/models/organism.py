from xml.etree import ElementTree

from django.db import models
from django.utils import timezone

import requests
from computedfields.models import ComputedFieldsModel, computed

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.base_models import TimestampedModel
from data_refinery_common.models.compendium_result import CompendiumResult
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


class UnknownOrganismId(Exception):
    pass


def get_scientific_name(taxonomy_id: int) -> str:
    parameters = {"db": TAXONOMY_DATABASE, "id": str(taxonomy_id), "api_key": NCBI_API_KEY}
    response = requests.get(EFETCH_URL, parameters)

    try:
        root = ElementTree.fromstring(response.text)
        taxon_list = root.findall("Taxon")
    except Exception:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(taxon_list) == 0:
        logger.error(
            "No names returned by ncbi.nlm.nih.gov for organism " + "with taxonomy ID %d.",
            taxonomy_id,
        )
        raise InvalidNCBITaxonomyId

    return taxon_list[0].find("ScientificName").text.upper()


def get_taxonomy_id(organism_name: str) -> int:
    parameters = {"db": TAXONOMY_DATABASE, "term": organism_name, "api_key": NCBI_API_KEY}
    response = requests.get(ESEARCH_URL, parameters)

    try:
        root = ElementTree.fromstring(response.text)
        id_list = root.find("IdList").findall("Id")
    except Exception:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(id_list) == 0:
        # If we can't find the taxonomy id, it's likely a bad organism name.
        error_message = (
            "Unable to retrieve NCBI taxonomy ID number for organism "
            + "with name: {}".format(organism_name)
        )
        logger.error(error_message)
        raise UnknownOrganismId(error_message)
    elif len(id_list) > 1:
        logger.warn(
            "Organism with name %s returned multiple NCBI taxonomy ID " + "numbers.", organism_name
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
    except Exception:
        logger.error("Bad response from eUtils.", text=response.text)
        raise

    if len(id_list) == 0:
        raise UnscientificNameError
    elif len(id_list) > 1:
        logger.warn(
            "Organism with name %s returned multiple NCBI taxonomy ID " + "numbers.", organism_name
        )

    return int(id_list[0].text)


class Organism(ComputedFieldsModel, TimestampedModel):
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

    experiments = models.ManyToManyField("Experiment", through="ExperimentOrganismAssociation")
    qn_target = models.ForeignKey(
        "ComputationalResult", blank=True, null=True, on_delete=models.SET_NULL
    )

    @computed(
        models.BooleanField(blank=False, null=True),
        depends=[["primary_compendium_results", ["quant_sf_only"]]],
    )
    def has_compendia(self):
        results = CompendiumResult.objects.all().filter(primary_organism=self, quant_sf_only=False)
        return results.count() != 0

    @computed(
        models.BooleanField(blank=False, null=True),
        depends=[["primary_compendium_results", ["quant_sf_only"]]],
    )
    def has_quantfile_compendia(self):
        results = CompendiumResult.objects.all().filter(primary_organism=self, quant_sf_only=True)
        return results.count() != 0

    def __str__(self):
        return str(self.name)

    def get_genus(self):
        return self.name.split("_")[0]

    @classmethod
    def get_or_create_object_for_id(cls, taxonomy_id: int):
        organism = (
            cls.objects.filter(taxonomy_id=taxonomy_id).order_by("-is_scientific_name").first()
        )

        if not organism:
            name = get_scientific_name(taxonomy_id).upper().replace(" ", "_")
            organism = Organism(name=name, taxonomy_id=taxonomy_id, is_scientific_name=True)
            organism.save()

        return organism

    @classmethod
    def get_name_for_id(cls, taxonomy_id: int) -> str:
        organism = cls.get_or_create_object_for_id(taxonomy_id)

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
                name=name, taxonomy_id=taxonomy_id, is_scientific_name=is_scientific_name
            )
            organism.save()

        return organism.taxonomy_id

    @classmethod
    def get_object_for_name(cls, name: str, taxonomy_id=None):
        name = name.upper()
        name = name.replace(" ", "_")
        organism = cls.objects.filter(name=name).first()

        if not organism:
            is_scientific_name = False
            if not taxonomy_id:
                try:
                    taxonomy_id = get_taxonomy_id_scientific(name)
                    is_scientific_name = True
                except UnscientificNameError:
                    taxonomy_id = get_taxonomy_id(name)

            organism = Organism(
                name=name, taxonomy_id=taxonomy_id, is_scientific_name=is_scientific_name
            )
            organism.save()

        return organism

    @classmethod
    def get_objects_with_qn_targets(cls):
        """ Return a list of Organisms who already have valid QN targets associated with them. """
        return Organism.objects.all().filter(qn_target__isnull=False)

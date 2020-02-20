import csv
import re
import shutil
import urllib
from abc import ABC
from typing import Dict, List

from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Organism, OriginalFile, SurveyJobKeyValue
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor

logger = get_and_configure_logger(__name__)


CHUNK_SIZE = 1024 * 256  # chunk_size is in bytes


MAIN_DIVISION_URL_TEMPLATE = "https://rest.ensembl.org/info/species?content-type=application/json"
DIVISION_URL_TEMPLATE = (
    "https://rest.ensembl.org/info/genomes/division/{division}?content-type=application/json"
)

SPECIES_DETAIL_URL_TEMPLATE = (
    "ftp://ftp.ensemblgenomes.org/pub/{short_division}/current/species_{division}.txt"
)
TRANSCRIPTOME_URL_TEMPLATE = (
    "ftp://ftp.{url_root}/fasta/{collection}{species_sub_dir}/dna/"
    "{filename_species}.{assembly}.dna.{schema_type}.fa.gz"
)
GTF_URL_TEMPLATE = (
    "ftp://ftp.{url_root}/gtf/{collection}{species_sub_dir}/"
    "{filename_species}.{assembly}.{assembly_version}.gtf.gz"
)


# For whatever reason the division in the download URL is shortened in
# a way that doesn't seem to be discoverable programmatically. I've
# therefore created this lookup map:
DIVISION_LOOKUP = {
    "EnsemblPlants": "plants",
    "EnsemblFungi": "fungi",
    "EnsemblBacteria": "bacteria",
    "EnsemblProtists": "protists",
    "EnsemblMetazoa": "metazoa",
}


# Ensembl will periodically release updated versions of the
# assemblies.  All divisions other than the main one have identical
# release versions. These urls will return what the most recent
# release version is.
MAIN_RELEASE_URL = "https://rest.ensembl.org/info/software?content-type=application/json"
DIVISION_RELEASE_URL = "https://rest.ensembl.org/info/eg_version?content-type=application/json"


def get_strain_mapping_for_organism(
    species_name: str, config_file="config/organism_strain_mapping.csv"
) -> List[Dict]:
    """Returns the row of the strain/organism mapping for the species_name
    """
    upper_name = species_name.upper()
    with open(config_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["organism"] == upper_name:
                return row

    return None


def get_species_detail_by_assembly(assembly: str, division: str) -> str:
    """Returns additional detail about a species given an assembly and a division.

    These details are necessary because the FTP directory for
    EnsemblBacteria and EnsemblFungi have an additional level in their
    paths that can only be determined by parsing this file. I found
    this out via the Ensembl dev mailing list.
    """
    bacteria_species_detail_url = SPECIES_DETAIL_URL_TEMPLATE.format(
        short_division=DIVISION_LOOKUP[division], division=division
    )

    urllib.request.urlcleanup()
    collection_path = assembly + "_collection.tsv"
    with open(collection_path, "wb") as collection_file:
        with urllib.request.urlopen(bacteria_species_detail_url) as request:
            shutil.copyfileobj(request, collection_file, CHUNK_SIZE)

    # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
    urllib.request.urlcleanup()

    with open(collection_path) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            if row["assembly"] == assembly:
                return row


class EnsemblUrlBuilder(ABC):
    """Generates URLs for different divisions of Ensembl.

    Each division of Ensembl has different conventions for its
    URLs. The logic contained in the init method of this base class is
    appropriate for most, but not all of the divisions. However, the
    logic contained in the build_* methods of this class is
    appropriate for all divisions.
    """

    def __init__(self, species: Dict):
        """Species is a Dict containing parsed JSON from the Division API."""
        self.url_root = "ensemblgenomes.org/pub/release-{assembly_version}/{short_division}"
        self.division = species["division"]
        self.short_division = DIVISION_LOOKUP[species["division"]]

        mapping = get_strain_mapping_for_organism(species["name"])
        if mapping:
            self.assembly = mapping["assembly"]
            self.strain = mapping["strain"]
        else:
            self.assembly = species["assembly_name"].replace(" ", "_")
            self.strain = None

        assembly_response = utils.requests_retry_session().get(DIVISION_RELEASE_URL)
        self.assembly_version = assembly_response.json()["version"]
        self.species_sub_dir = species["name"]
        self.filename_species = species["name"].capitalize()

        # These fields aren't needed for the URL, but they vary between
        # the two REST APIs.
        self.scientific_name = species["name"].upper()

        # This field can be stored in multiple keys, but if
        # `species_taxonomy_id` is there it's the one we want because
        # it's not strain-specific.
        if "species_taxonomy_id" in species:
            self.taxonomy_id = species["species_taxonomy_id"]
        else:
            self.taxonomy_id = species["taxonomy_id"]

        # This field is only needed for EnsemblBacteria and EnsemblFungi.
        self.collection = ""

    def build_transcriptome_url(self) -> str:
        url_root = self.url_root.format(
            assembly_version=self.assembly_version, short_division=self.short_division
        )
        url = TRANSCRIPTOME_URL_TEMPLATE.format(
            url_root=url_root,
            species_sub_dir=self.species_sub_dir,
            collection=self.collection,
            filename_species=self.filename_species,
            assembly=self.assembly,
            schema_type="primary_assembly",
        )

        # If the primary_assembly is not available use toplevel instead.
        try:
            # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
            urllib.request.urlcleanup()
            file_handle = urllib.request.urlopen(url)
            file_handle.close()
            urllib.request.urlcleanup()
        except Exception:
            url = url.replace("primary_assembly", "toplevel")

        return url

    def build_gtf_url(self) -> str:
        url_root = self.url_root.format(
            assembly_version=self.assembly_version, short_division=self.short_division
        )
        return GTF_URL_TEMPLATE.format(
            url_root=url_root,
            species_sub_dir=self.species_sub_dir,
            collection=self.collection,
            filename_species=self.filename_species,
            assembly=self.assembly,
            assembly_version=self.assembly_version,
        )


class MainEnsemblUrlBuilder(EnsemblUrlBuilder):
    """Special logic specific to the main Ensembl division.

    There is one Ensembl division which is just called Ensembl. This
    is confusing so I refer to it as the main Ensembl division. It
    follows the same general pattern as the rest of them for URLs, but
    just not quite the same base URL structure. Also its REST API
    returns JSON with similar data except with slightly different key
    names.
    """

    def __init__(self, species: Dict):
        self.url_root = "ensembl.org/pub/release-{assembly_version}"
        self.short_division = None
        self.species_sub_dir = species["name"]
        self.collection = ""
        self.filename_species = species["name"].capitalize()
        self.assembly = species["assembly"]
        self.assembly_version = (
            utils.requests_retry_session().get(MAIN_RELEASE_URL).json()["release"]
        )
        self.scientific_name = self.filename_species.replace("_", " ")
        self.taxonomy_id = species["taxon_id"]


class EnsemblProtistsUrlBuilder(EnsemblUrlBuilder):
    """Special logic specific to the EnsemblProtists division.

    EnsemblProtists is special because the first letter of the species
    name is always capitalized within the name of the file, instead of
    only when there's not a collection subnested.
    """

    def __init__(self, species: Dict):
        super().__init__(species)
        self.filename_species = species["name"].capitalize()


class EnsemblBacteriaUrlBuilder(EnsemblUrlBuilder):
    """The EnsemblBacteria URLs are extra tricky because they have an extra layer in them.

    This requires parsing a special file to find out what collection
    the species belongs to.
    """

    def __init__(self, species: Dict):
        super().__init__(species)
        species_detail = get_species_detail_by_assembly(self.assembly, self.division)

        if species_detail:
            self.species_sub_dir = species_detail["species"]
            self.filename_species = species_detail["species"].capitalize()
            collection_pattern = r"bacteria_.+_collection"
            match = re.match(collection_pattern, species_detail["core_db"])
            if match:
                # Need to append a / to the collection because it's not
                # present in all the routes so we don't want to put it in the
                # template and end up with a // in the path if collection is
                # blank.
                self.collection = match.group(0) + "/"


class EnsemblFungiUrlBuilder(EnsemblUrlBuilder):
    """The EnsemblFungi URLs work similarly to EnsemblBacteria.

    EnsemblFungi is special because there is an assembly_name TIGR
    which needs to be corrected to CADRE for some reason.
    """

    def __init__(self, species: Dict):
        super().__init__(species)
        if self.assembly == "TIGR":
            self.assembly = "CADRE"

        species_detail = get_species_detail_by_assembly(self.assembly, self.division)

        if species_detail:
            self.species_sub_dir = species_detail["species"]
            self.filename_species = species_detail["species"].capitalize()
            collection_pattern = r"fungi_.+_collection"
            match = re.match(collection_pattern, species_detail["core_db"])
            if match:
                # Need to append a / to the collection because it's not
                # present in all the routes so we don't want to put it in the
                # template and end up with a // in the path if collection is
                # blank.
                self.collection = match.group(0) + "/"


def ensembl_url_builder_factory(species: Dict) -> EnsemblUrlBuilder:
    """Returns instance of EnsemblUrlBuilder or one of its subclasses.

    The class of the returned object is based on the species' division.
    """
    if species["division"] == "EnsemblProtists":
        return EnsemblProtistsUrlBuilder(species)
    elif species["division"] == "EnsemblFungi":
        return EnsemblFungiUrlBuilder(species)
    elif species["division"] == "EnsemblVertebrates":
        return MainEnsemblUrlBuilder(species)
    elif species["division"] == "EnsemblBacteria":
        return EnsemblBacteriaUrlBuilder(species)
    else:
        return EnsemblUrlBuilder(species)


class TranscriptomeIndexSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.TRANSCRIPTOME_INDEX.value

    def _generate_files(self, species: Dict) -> None:
        url_builder = ensembl_url_builder_factory(species)
        fasta_download_url = url_builder.build_transcriptome_url()
        gtf_download_url = url_builder.build_gtf_url()

        # Getting the object will ensure it is created in the DB.
        Organism.get_or_create_object_for_id(url_builder.taxonomy_id)

        all_new_files = []

        fasta_filename = url_builder.filename_species + ".fa.gz"
        original_file = OriginalFile()
        original_file.source_filename = fasta_filename
        original_file.source_url = fasta_download_url
        original_file.is_archive = True
        original_file.is_downloaded = False
        original_file.save()
        all_new_files.append(original_file)

        gtf_filename = url_builder.filename_species + ".gtf.gz"
        original_file = OriginalFile()
        original_file.source_filename = gtf_filename
        original_file.source_url = gtf_download_url
        original_file.is_archive = True
        original_file.is_downloaded = False
        original_file.save()
        all_new_files.append(original_file)

        return all_new_files

    def survey(self, source_type=None) -> bool:
        """
        Surveying here is a bit different than discovering an experiment
        and samples.
        """
        if source_type != "TRANSCRIPTOME_INDEX":
            return False

        try:
            species_files = self.discover_species()
        except Exception:
            logger.exception(
                "Exception caught while discovering species. Terminating survey job.",
                survey_job=self.survey_job.id,
            )
            return False

        try:
            for specie_file_list in species_files:
                self.queue_downloader_job_for_original_files(
                    specie_file_list, is_transcriptome=True
                )
        except Exception:
            logger.exception(
                "Failed to queue downloader jobs. Terminating survey job.",
                survey_job=self.survey_job.id,
            )
            return False

        return True

    def discover_species(self):
        ensembl_division = SurveyJobKeyValue.objects.get(
            survey_job_id=self.survey_job.id, key__exact="ensembl_division"
        ).value

        logger.info(
            "Surveying %s division of ensembl.", ensembl_division, survey_job=self.survey_job.id,
        )

        try:
            organism_name = SurveyJobKeyValue.objects.get(
                survey_job_id=self.survey_job.id, key__exact="organism_name"
            ).value
            organism_name = organism_name.lower().replace(" ", "_")
        except SurveyJobKeyValue.DoesNotExist:
            organism_name = None

        if ensembl_division in ["EnsemblFungi", "EnsemblBacteria"]:
            if organism_name is None:
                logger.error(
                    "Organism name must be specified for Fungi and Bacteria divisions.",
                    ensembl_division=ensembl_division,
                    organism_name=organism_name,
                )
                return []
            else:
                if get_strain_mapping_for_organism(organism_name) is None:
                    logger.error(
                        (
                            "Organism name must be listed in config/organism_strain_"
                            "mappings.csv for Fungi and Bacteria divisions."
                        ),
                        ensembl_division=ensembl_division,
                        organism_name=organism_name,
                    )
                    return []

        # The main division has a different base URL for its REST API.
        if ensembl_division == "Ensembl":
            r = utils.requests_retry_session().get(MAIN_DIVISION_URL_TEMPLATE)

            # Yes I'm aware that specieses isn't a word. However I need to
            # distinguish between a singlular species and multiple species.
            specieses = r.json()["species"]
        else:
            formatted_division_url = DIVISION_URL_TEMPLATE.format(division=ensembl_division)
            r = utils.requests_retry_session().get(formatted_division_url)
            specieses = r.json()

        all_new_species = []
        if organism_name:
            for species in specieses:
                if ensembl_division == "EnsemblFungi" and organism_name in species["name"]:
                    # Fungi have a strain identifier in their
                    # names. This is different than everything else,
                    # so we're going to handle this special case by
                    # just overwriting this. This is okay because we
                    # just have to discover one species for the
                    # organism, and then our strain mapping will make
                    # sure we use the correct strain and assembly.
                    species["name"] = organism_name

                    all_new_species.append(self._generate_files(species))
                    break
                elif "name" in species and organism_name == species["name"]:
                    all_new_species.append(self._generate_files(species))
                    break
        else:
            for species in specieses:
                all_new_species.append(self._generate_files(species))

        if len(all_new_species) == 0:
            logger.error(
                "Unable to find any species!",
                ensembl_division=ensembl_division,
                organism_name=organism_name,
            )

        return all_new_species

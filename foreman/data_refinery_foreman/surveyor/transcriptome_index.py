import re
import urllib
from abc import ABC
from typing import Dict, List

from django.utils import timezone

from data_refinery_common.job_lookup import Downloaders, ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OriginalFile, SurveyJobKeyValue
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor

logger = get_and_configure_logger(__name__)


MAIN_DIVISION_URL_TEMPLATE = "https://rest.ensembl.org/info/species?content-type=application/json"
DIVISION_URL_TEMPLATE = (
    "https://rest.ensembl.org/info/genomes/division/{division}" "?content-type=application/json"
)

TRANSCRIPTOME_URL_TEMPLATE = (
    "ftp://ftp.{url_root}/fasta/{species_sub_dir}/dna/"
    "{filename_species}.{assembly}.dna.{schema_type}.fa.gz"
)
GTF_URL_TEMPLATE = (
    "ftp://ftp.{url_root}/gtf/{species_sub_dir}/"
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
        self.short_division = DIVISION_LOOKUP[species["division"]]
        self.assembly = species["assembly_name"].replace(" ", "_")
        self.assembly_version = (
            utils.requests_retry_session().get(DIVISION_RELEASE_URL).json()["version"]
        )

        self.species_sub_dir = species["name"]
        self.filename_species = species["name"].capitalize()

        # These fields aren't needed for the URL, but they vary between
        # the two REST APIs.
        self.scientific_name = species["name"].upper()
        self.taxonomy_id = species["taxonomy_id"]

    def build_transcriptome_url(self) -> str:
        url_root = self.url_root.format(
            assembly_version=self.assembly_version, short_division=self.short_division
        )
        url = TRANSCRIPTOME_URL_TEMPLATE.format(
            url_root=url_root,
            species_sub_dir=self.species_sub_dir,
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
        except:
            url = url.replace("primary_assembly", "toplevel")

        return url

    def build_gtf_url(self) -> str:
        url_root = self.url_root.format(
            assembly_version=self.assembly_version, short_division=self.short_division
        )
        return GTF_URL_TEMPLATE.format(
            url_root=url_root,
            species_sub_dir=self.species_sub_dir,
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
        self.filename_species = species["species"].capitalize()


class EnsemblFungiUrlBuilder(EnsemblProtistsUrlBuilder):
    """The EnsemblFungi URLs work the similarly to Protists division.

    EnsemblFungi is special because there is an assembly_name TIGR
    which needs to be corrected to CADRE for some reason.
    """

    def __init__(self, species: Dict):
        super().__init__(species)
        if self.assembly == "TIGR":
            self.assembly = "CADRE"


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
    else:
        return EnsemblUrlBuilder(species)


class TranscriptomeIndexSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.TRANSCRIPTOME_INDEX.value

    def _clean_metadata(self, species: Dict) -> Dict:
        """Removes fields from metadata which shouldn't be stored.

        Also cast any None values to str so they can be stored in the
        database.
        These fields shouldn't be stored because:
        The taxonomy id is stored as fields on the Organism.
        Aliases and groups are lists we don't need.
        """
        species.pop("taxon_id") if "taxon_id" in species else None
        species.pop("taxonomy_id") if "taxonomy_id" in species else None
        species.pop("aliases") if "aliases" in species else None
        species.pop("groups") if "groups" in species else None

        # Cast to List since we're modifying the size of the dict
        # while iterating over it
        for k, v in list(species.items()):
            if v is None:
                species.pop(k)
            else:
                species[k] = str(v)

        return species

    def _generate_files(self, species: Dict) -> None:
        url_builder = ensembl_url_builder_factory(species)
        fasta_download_url = url_builder.build_transcriptome_url()
        gtf_download_url = url_builder.build_gtf_url()

        platform_accession_code = species.pop("division")
        self._clean_metadata(species)

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
                ("Exception caught while discovering species. " "Terminating survey job."),
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
                ("Failed to queue downloader jobs. " "Terminating survey job."),
                survey_job=self.survey_job.id,
            )
            return False

        return True

    def discover_species(self):
        ensembl_division = SurveyJobKeyValue.objects.get(
            survey_job_id=self.survey_job.id, key__exact="ensembl_division"
        ).value

        logger.info(
            "Surveying %s division of ensembl.", ensembl_division, survey_job=self.survey_job.id
        )

        # The main division has a different base URL for its REST API.
        if ensembl_division == "Ensembl":
            r = utils.requests_retry_session().get(MAIN_DIVISION_URL_TEMPLATE)

            # Yes I'm aware that specieses isn't a word. However I need to
            # distinguish between a singlular species and multiple species.
            specieses = r.json()["species"]
        else:
            r = utils.requests_retry_session().get(
                DIVISION_URL_TEMPLATE.format(division=ensembl_division)
            )
            specieses = r.json()

        try:
            organism_name = SurveyJobKeyValue.objects.get(
                survey_job_id=self.survey_job.id, key__exact="organism_name"
            ).value
            organism_name = organism_name.lower().replace(" ", "_")
        except SurveyJobKeyValue.DoesNotExist:
            organism_name = None

        all_new_species = []
        if organism_name:
            for species in specieses:
                # This key varies based on whether the division is the
                # main one or not... why couldn't they just make them
                # consistent?
                if ("species" in species and species["species"] == organism_name) or (
                    "name" in species and species["name"] == organism_name
                ):
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

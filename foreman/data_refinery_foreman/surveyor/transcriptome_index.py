from abc import ABC
import requests
import re
from pprint import pprint
import urllib
from typing import List, Dict
from django.utils import timezone
from data_refinery_common.models import (
    Batch,
    File,
    SurveyJobKeyValue
)
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


DIVISION_URL_TEMPLATE = ("http://rest.ensemblgenomes.org/info/genomes/division/{division}"
                         "?content-type=application/json")
TRANSCRIPTOME_URL_TEMPLATE = (
    "ftp://ftp.ensemblgenomes.org/pub/release-37/{short_division}/fasta/"
    "{species_sub_dir}/dna/{file_name_species}.{assembly}.dna.{schema_type}.fa.gz"
)
GTF_URL_TEMPLATE = ("ftp://ftp.ensemblgenomes.org/pub/release-37/{short_division}/gtf/"
                    "{species_sub_dir}/{file_name_species}.{assembly}.37.gtf.gz")


# For whatever reason the division in the download URL is shortened in
# a way that doesn't seem to be discoverable programmatically. I've
# therefore created this lookup map:
DIVISION_LOOKUP = {"EnsemblPlants": "plants",
                   "EnsemblFungi": "fungi",
                   "EnsemblBacteria": "bacteria",
                   "EnsemblProtists": "protists",
                   "EnsemblMetazoa": "metazoa"}


class EnsemblUrlLogic(ABC):
    """Each division of Ensembl has different conventions for its URLs.

    The logic contained in the class is appropriate for most, but not
    all of them. All of the logic is performed in the init method, and
    then the results can be found by accessing the following member
    variables:
    short_division: The division name as it appears in URLs.
    assembly: The name of the assembly of the genome?
    species_sub_dir: The correct sub-directory for the species.
    file_name_species: The species' name as it appears in the file name.
    """

    def __init__(self, species: Dict):
        """Species is a Dict containing parsed JSON from the Division API."""
        self.short_division = DIVISION_LOOKUP[species["division"]]
        self.assembly = species["assembly_name"].replace(" ", "_")

        # Some species are nested within a collection directory. If
        # this is the case, then we need to add that extra directory
        # to the URL, and for whatever reason the filename is not
        # capitalized.
        COLLECTION_REGEX = r"^(.*_collection).*"
        match_object = re.search(COLLECTION_REGEX, species["dbname"])
        if match_object:
            self.species_sub_dir = match_object.group(1) + "/" + species["species"]
            self.file_name_species = species["species"]
        else:
            self.species_sub_dir = species["species"]
            self.file_name_species = species["species"].capitalize()

    # def get_species_sub_dir(self) -> str:
    #     """Returns the correct sub-directory for the species."""
    #     return

    # def get_file_name_species(self) -> str:
    #     """Returns the species' name as it appears in the file name."""
    #     return

    # def get_short_division(self) -> str:
    #     """Returns a shorter division name.

    #     For whatever reason the division in the download URL is
    #     shortened in a way that doesn't seem to be discoverable
    #     programmatically so I've hardcoded them like this.
    #     """
    #     return


class EnsemblProtistsUrlLogic(EnsemblUrlLogic):
    """Special logic specific to the EnsemblProtists division.

    EnsemblProtists is special because the first letter of the species
    name is always capitalized within the name of the file, instead of
    only when there's not a collection subnested.
    """

    def __init__(self, species: Dict):
        super().__init__(species)
        self.file_name_species = species["species"].capitalize()

    # def get_species_sub_dir(self) -> str:
    #     return self.species["species"]

    # def get_file_name_species(self) -> str:
    #     return self.species["species"].capitalize()

    # def get_short_division(self) -> str:
    #     return "plants"


# class EnsemblFungiUrlLogic(EnsemblUrlLogic):
#     """Special logic specific to the EnsemblFungi division.

#     EnsemblPlants is special because the first letter of the
#     species name gets capitalized within the name of the file.
#     """

#     def get_species_sub_dir(self) -> str:
#         COLLECTION_REGEX = r"^(.*_collection).*"
#         match_object = re.search(COLLECTION_REGEX, self.species["dbname"])
#         if match_object:
#             return match_object.group(1) + "/" + self.species["species"]
#         else:
#             return self.species["species"]

#     def get_file_name_species(self) -> str:
#         return self.species["species"].capitalize()

#     def get_short_division(self) -> str:
#         return "plants"


def ensembl_url_logic_factory(species: Dict) -> EnsemblUrlLogic:
    if species["division"] == "EnsemblProtists":
        return EnsemblProtistsUrlLogic(species)
    # elif species["division"] == "EnsemblFungi":
    #     return EnsemblFungiUrlLogic(species)
    # elif species["division"] == "EnsemblBacteria":
    #     return EnsemblBacteriaUrlLogic(species)
    # elif species["division"] == "EnsemblProtists":
    #     return EnsemblProtistsUrlLogic(species)
    # elif species["division"] == "EnsemblMetazoa":
    #     return EnsemblMetazoaUrlLogic(species)
    else:
        return EnsemblUrlLogic(species)


class TranscriptomeIndexSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.TRANSCRIPTOME_INDEX.value

    def determine_pipeline(self, batch: Batch, key_values: Dict = {}) -> ProcessorPipeline:
        return ProcessorPipeline.TRANSCRIPTOME_INDEX

    def group_batches(self) -> List[List[Batch]]:
        """Groups batches based on the download URL of their only File.

        This is a duplicate of array_express' group_batches fn. Can I
        DRY this up?
        """
        # Builds a mapping of each unique download_url to a list of
        # Batches whose first File's download_url matches.
        download_url_mapping = {}
        for batch in self.batches:
            download_url = batch.files[0].download_url
            if download_url in download_url_mapping:
                download_url_mapping[download_url].append(batch)
            else:
                download_url_mapping[download_url] = [batch]

        # The values of the mapping we built are the groups the
        # batches should be grouped into.
        return list(download_url_mapping.values())

    # def _plant_special_logic(self, species: Dict) -> Tuple:
    #     """Special logic specific to the EnsemblPlants division.

    #     EnsemblPlants is special because the first letter of the
    #     species name gets capitalized within the name of the file.
    #     """

    def _generate_batches(self, species: Dict) -> None:
        pprint(species)

        # COLLECTION_REGEX = r"^(.*_collection).*"
        # match_object = re.search(COLLECTION_REGEX, species["dbname"])
        # if match_object:
        #     species_sub_dir = match_object.group(1) + "/" + species["species"]
        #     file_name_species = species["species"]
        # else:
        #     species_sub_dir = species["species"]
        #     file_name_species = species["species"].capitalize()

        url_fields = ensembl_url_logic_factory(species)

        # assembly = species["assembly_name"].replace(" ", "_")

        # short_division = DIVISION_LOOKUP[species["division"]]

        fasta_download_url = TRANSCRIPTOME_URL_TEMPLATE.format(
            short_division=url_fields.short_division,
            species_sub_dir=url_fields.species_sub_dir,
            file_name_species=url_fields.file_name_species,
            assembly=url_fields.assembly,
            schema_type="primary_assembly")
        # If the primary_assembly is not available use toplevel instead.
        try:
            file_handle = urllib.request.urlopen(fasta_download_url)
            file_handle.close()
        except:
            fasta_download_url = fasta_download_url.replace("primary_assembly", "toplevel")

        gtf_download_url = GTF_URL_TEMPLATE.format(
            short_division=url_fields.short_division,
            species_sub_dir=url_fields.species_sub_dir,
            file_name_species=url_fields.file_name_species,
            assembly=url_fields.assembly)

        current_time = timezone.now()
        for length in ("_long", "_short"):
            fasta_file = File(name=species["species"] + length + ".fa.gz",
                              download_url=fasta_download_url,
                              raw_format="fa.gz",
                              processed_format="tar.gz",
                              size_in_bytes=-1)  # Will have to be determined later

            gtf_file = File(name=species["species"] + length + ".gtf.gz",
                            download_url=gtf_download_url,
                            raw_format="gtf.gz",
                            processed_format="tar.gz",
                            size_in_bytes=-1)  # Will have to be determined later

            self.add_batch(platform_accession_code=species["division"],
                           experiment_accession_code=species["species"],
                           organism_id=species["taxonomy_id"],
                           organism_name=species["name"].upper(),
                           experiment_title="NA",
                           release_date=current_time,
                           last_uploaded_date=current_time,
                           files=[fasta_file, gtf_file],
                           # Store the rest of the metadata about these!
                           key_values={"length": length})

    def discover_batches(self):
        ensembl_division = (
            SurveyJobKeyValue
            .objects
            .get(survey_job_id=self.survey_job.id,
                 key__exact="ensembl_division")
            .value
        )

        logger.info("Surveying %s division of ensembl.",
                    ensembl_division,
                    survey_job=self.survey_job.id)

        r = requests.get(DIVISION_URL_TEMPLATE.format(division=ensembl_division))
        # Yes I'm aware that specieses isn't a word. However I need to
        # distinguish between a singlular species and multiple species.
        specieses = r.json()
        # pprint(specieses)

        for species in specieses:
            self._generate_batches(species)
            # TEMPORARY FOR TESTING
            break

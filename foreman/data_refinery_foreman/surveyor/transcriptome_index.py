import requests
from typing import List, Dict
from django.utils import timezone
from data_refinery_common.models import (
    Batch,
    File,
    SurveyJobKeyValue,
    Organism
)
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


DIVISION_URL_TEMPLATE = ("http://rest.ensemblgenomes.org/info/genomes/division/{division}"
                         "?content-type=application/json")
TRANSCRIPTOME_URL_TEMPLATE = (
    "ftp://ftp.ensemblgenomes.org/pub/release-37/{short_division}/fasta/"
    "{species}/dna/{caps_species}.{assembly}.dna.toplevel.fa.gz"
)
GTF_URL_TEMPLATE = ("ftp://ftp.ensemblgenomes.org/pub/release-37/{short_division}/gtf/"
                    "{species}/{caps_species}.{assembly}.37.gtf.gz")


# For whatever reason the division in the download URL is shortened in
# a way that doesn't seem to be discoverable programmatically. I've
# therefore created this lookup map:
DIVISION_LOOKUP = {"EnsemblPlants": "plants",
                   "EnsemblFungi": "fungi",
                   "EnsemblBacteria": "bacteria",
                   "EnsemblProtists": "protists",
                   "EnsemblMetazoa": "metazoa"}


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

    def _generate_batches(self, species: Dict) -> None:
        # The caps_species is a string which is formatted like this
        # because that's how it is in the URL:
        caps_species = species["name"].replace(" ", "_")
        short_division = DIVISION_LOOKUP[species["division"]]
        fasta_download_url = TRANSCRIPTOME_URL_TEMPLATE.format(short_division=short_division,
                                                               species=species["species"],
                                                               caps_species=caps_species,
                                                               assembly=species["assembly_name"])
        gtf_download_url = GTF_URL_TEMPLATE.format(short_division=short_division,
                                                   species=species["species"],
                                                   caps_species=caps_species,
                                                   assembly=species["assembly_name"])

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

        ensembl_division = "EnsemblPlants"

        r = requests.get(DIVISION_URL_TEMPLATE.format(division=ensembl_division))
        # Yes I'm aware that specieses isn't a word. However I need to
        # distinguish between a singlular species and multiple species.
        specieses = r.json()

        for species in specieses:
            self._generate_batches(species)
            # TEMPORARY FOR TESTING
            break

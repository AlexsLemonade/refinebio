"""
RNA-Seq accession gathering automation.
Data source: https://www.ebi.ac.uk/ena/portal/api/
"""

from json.decoder import JSONDecodeError
from typing import List, Set
from urllib.parse import quote

import requests
from retrying import retry

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.gathered_accession import GatheredAccession
from data_refinery_foreman.gatherer.agents.base import AccessionAgentBase

logger = get_and_configure_logger(__name__)


class RNASeqAgent(AccessionAgentBase):
    """
    RNA-Seq accession gathering agent. The data is fetched from
    The European Nucleotide Archive (ENA) Portal.
    See https://www.ebi.ac.uk/ena/portal/api/ for more information about the API
    endpoints.
    """

    DATA_CHUNK_SIZE = 10000
    DATA_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
    SOURCE = "ebi-ena-portal"
    SOURCE_NAME = "rna-seq"
    TECHNOLOGY = "rna-seq"

    def __str__(self):
        return "RNA-Seq accession agent"

    def build_query(self, taxon_id: str = None) -> str:
        """
        Returns a query to use for getting specific taxon ID accessions.
        Some special characters must remain unquoted.
        """

        AND = " AND "
        OR = " OR "
        # TODO(ark): extract instrument models to a config file.
        instrument_models = (
            "HiSeq X Five",
            "HiSeq X Ten",
            "Illumina Genome Analyzer II",
            "Illumina Genome Analyzer IIx",
            "Illumina Genome Analyzer",
            "Illumina HiScanSQ",
            "Illumina HiSeq 1000",
            "Illumina HiSeq 1500",
            "Illumina HiSeq 2000",
            "Illumina HiSeq 2500",
            "Illumina HiSeq 3000",
            "Illumina HiSeq 4000",
            "Illumina MiSeq",
            "Illumina NovaSeq 6000",
            "Ion Torrent Proton",
            "Ion Torrent S5 XL",
            "Ion Torrent S5",
            "NextSeq 500",
            "NextSeq 550",
        )

        instrument_models = OR.join((f'instrument_model="{im}"' for im in instrument_models))
        conditions = [
            # Relevant date fields: collection_date, collection_date_submitted,
            # first_public, last_updated.
            f"first_public >= {self.since}",
            f"first_public <= {self.until}",
            f"({instrument_models})",
            'library_source="TRANSCRIPTOMIC"',
            'library_strategy="RNA-Seq"',
        ]

        if taxon_id:
            conditions.append(f"tax_eq({taxon_id})")
        elif self.keyword:
            search_fields = (
                "assembly_software",
                "bio_material",
                "center_name",
                "collected_by",
                "experiment_title",
                "host_body_site",
                "instrument_model",
                "instrument_platform",
                "library_name",
                "project_name",
                "sample_title",
                "sequencing_method",
                "study_title",
            )
            search_fields = OR.join(
                (f'{sf}="*{self.keyword}*"' for sf in search_fields)
            )  # Keyword regex.
            conditions.append(f"({search_fields})")
        elif self.organism:
            # `host`: Natural (as opposed to laboratory) host to the organism from which sample
            #         was obtained.
            # `host_scientific_name`: Scientific name of the natural (as opposed to laboratory)
            #                         host to the organism from which sample was obtained.
            # `scientific_name` Scientific name of the organism from which the sample was derived.
            # Neither `host_scientific_name` nor `scientific_name` available for search.
            # https://www.ebi.ac.uk/ena/portal/api/searchFields?dataPortal=ena&format=json&result=read_study
            conditions.append(f'host="{self.organism}"')

        return quote(AND.join(conditions), safe='*()-="<>/ ')  # Must remain unquoted.

    def collect_data(self) -> Set[str]:
        """Gets new accessions from EBI ENA API."""
        accessions = set()

        if self.ids:
            logger.debug(
                f"Getting RNA-Seq entries by taxon ID(s): "
                f"{', '.join((str(i) for i in self.ids))} for [{self.since} - {self.until}] range."
            )
            total = len(self.ids)
            for idx, taxon_id in enumerate(self.ids):
                if self.count and len(accessions) >= self.count:
                    break

                if total > 1:
                    logger.debug(f"Getting entries for taxon ID {taxon_id}, {idx + 1} of {total}.")
                accessions.update(self.fetch_data(taxon_id=taxon_id))
        elif self.keyword:
            logger.debug(
                f'Getting RNA-Seq entries by "{self.keyword}" keyword '
                f"for [{self.since} - {self.until}] range."
            )
            accessions.update(self.fetch_data())
        elif self.organism:
            logger.debug(
                f'Getting entries by "{self.organism}" organism '
                f"for [{self.since} - {self.until}] range."
            )
            accessions.update(self.fetch_data())

        return accessions

    def fetch_data(self, taxon_id=None) -> Set[str]:
        """
        Retrieves accessions from API search endpoint.
        The API allows to set limit to 0 (get all in one request) but we do
        it in a paginated fashion with `self.DATA_CHUNK_SIZE` as a page size.
        """

        @retry(**self.retry_params)
        def get_response(url, **kwargs):
            """Gets response from an API endpoint."""
            return requests.post(url, **kwargs)

        accessions = set()

        fields = [
            "first_public",
            "scientific_name",
            "secondary_study_accession",
        ]  # For DRP/ERP/SRP-prefixed accessions.
        data = {
            "dataPortal": "ena",
            # TODO(ark): add excludeAccessions/excludeAccessionType support.
            "fields": ",".join(fields),  # Use "all" to get all fields.
            "format": "json",
            "limit": self.DATA_CHUNK_SIZE,
            "offset": 0,
            "query": self.build_query(taxon_id=taxon_id),
            "result": "read_study",
            "sortFields": fields,
        }

        is_done = False
        while not is_done:
            logger.debug(
                f"Processing entries {data['offset'] + 1} - {data['offset'] + self.DATA_CHUNK_SIZE}"
            )
            entries = ()
            try:
                response = get_response(self.DATA_URL, data=data)
                entries = response.json()
                entries = (
                    GatheredAccession.create_from_external_entry(
                        entry, self.SOURCE, self.TECHNOLOGY
                    )
                    for entry in entries
                )
            except JSONDecodeError:
                is_done = True
            except TypeError:
                logger.error(f"Couldn't get data from {self.data_url}. Response: {entries}")
            data["offset"] += self.DATA_CHUNK_SIZE

            if self.previous_accessions:
                entries = (
                    entry
                    for entry in entries
                    if entry.accession_code not in self.previous_accessions
                )
            accessions.update(entries)

            # Quit after getting a sufficient amount of accessions.
            if self.count and len(accessions) >= self.count:
                is_done = True

        return accessions

    def get_ids(self) -> List[str]:
        """Returns a combined list of passed taxon IDs."""
        ids = set()

        if self.options["taxon_id"]:
            ids.update(self.options["taxon_id"])

        if self.options["taxon_ids_file"]:
            with open(self.options["taxon_ids_file"]) as taxon_id_file:
                ids.update((taxon_id.strip() for taxon_id in taxon_id_file.readlines()))

        return sorted(ids)

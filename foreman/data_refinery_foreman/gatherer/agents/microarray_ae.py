"""MicroArray ArrayExpress accession gathering automation.
Data source: https://www.ebi.ac.uk/biostudies/help"""

from typing import List, Set

import requests
from retrying import retry

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.gathered_accession import GatheredAccession
from data_refinery_foreman.gatherer.agents.base import AccessionAgentBase

logger = get_and_configure_logger(__name__)


class MicroArrayExpressAccessionAgent(AccessionAgentBase):
    """
    MicroArray ArrayExpress accession gathering agent. The data is fetched from
    the BioStudies database. See https://www.ebi.ac.uk/biostudies/help and
    https://www.ebi.ac.uk/biostudies/arrayexpress/help#programmatic for more
    information about the API endpoints.
    """

    DATA_CHUNK_SIZE = 100
    DATA_URL = "https://www.ebi.ac.uk/biostudies/api/v1/search"

    def __str__(self):
        return "MicroArray ArrayExpress accession agent"

    def build_query(self) -> dict:
        """Returns a query dict for getting array/organism specific accessions."""
        query_dict = {
            "directsub": "true",
            "page": 1,
            "pageSize": self.DATA_CHUNK_SIZE,
            "release_date": f"[{self.since} TO {self.until}]",
            "type": "study",
        }

        if self.ids:
            # TODO(ark): figure out better way of array filtering.
            # Also make sure it's equivalent to the array filtering in this query
            # https://github.com/AlexsLemonade/accession_retrieval/blob/master/experiment_accession_retrieval.R#L208
            query_dict.update({"content": ", ".join(self.ids)})
        elif self.keyword:
            query_dict.update({"content": self.keyword})
        elif self.organism:
            query_dict.update({"organism": f'"{self.organism}"'})

        return query_dict

    def collect_data(self) -> Set[str]:
        """Gets new accessions from EBI Biostudies API."""
        accessions = set()

        if self.ids:
            message = (
                "Getting MicroArray ArrayExpress entries by "
                f"ArrayExpress ID(s): {', '.join(self.ids)} for [{self.since} - {self.until}] "
                "range."
            )
        elif self.keyword:
            message = (
                "Getting MicroArray ArrayExpress entries by "
                f'"{self.keyword}" keyword for [{self.since} - {self.until}] range.'
            )
        elif self.organism:
            message = (
                "Getting MicroArray ArrayExpress entries by "
                f'"{self.organism}" organism for [{self.since} - {self.until}] range.'
            )
        else:
            return accessions

        logger.debug(message)
        accessions.update(self.fetch_data())

        return accessions

    def fetch_data(self) -> Set[str]:
        """Retrieves accessions from API search endpoint."""

        @retry(**self.retry_params)
        def get_response(url, **kwargs):
            """Gets response from an API endpoint."""
            return requests.get(url, **kwargs)

        accessions = set()

        is_done = False
        params = self.build_query()
        while not is_done:
            range_start = (params["page"] - 1) * params["pageSize"] + 1
            range_end = (params["page"] - 1) * params["pageSize"] + self.DATA_CHUNK_SIZE
            logger.debug(f"Processing entries {range_start} - {range_end}")

            response = get_response(self.DATA_URL, params=params)
            entries = response.json().get("hits", ())
            if entries:
                entries = (
                    GatheredAccession.create_from_ma_ae_entry(entry, organism=self.organism)
                    for entry in entries
                )
                params["page"] += 1
            else:
                is_done = True

            if self.previous_accessions:
                entries = (entry for entry in entries if entry.code not in self.previous_accessions)
            accessions.update(entries)

            # Quit after getting a sufficient amount of accessions.
            if self.count and len(accessions) >= self.count:
                is_done = True

        return accessions

    def get_ids(self) -> List[str]:
        """Returns a combined list of passed ArrayExpress IDs."""
        ids = set()

        if self.options["ae_id"]:
            ids.update(self.options["ae_id"])

        if self.options["ae_ids_file"]:
            with open(self.options["ae_ids_file"]) as ae_ids_file:
                ids.update((ae_id.strip() for ae_id in ae_ids_file.readlines()))

        return sorted(ids)

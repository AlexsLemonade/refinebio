"""
MicroArray GEO accession gathering automation.
Data source: local SQLite meta DB from
https://www.bioconductor.org/packages/release/bioc/html/GEOmetadb.html
"""

import os
import re
import sqlite3
from typing import List, Set

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.gathered_accession import GatheredAccession
from data_refinery_foreman.gatherer.agents.base import AccessionAgentBase

logger = get_and_configure_logger(__name__)


class GEOAgent(AccessionAgentBase):
    """
    MicroArray GEO accession gathering agent. The data is fetched from a local
    SQLite GEO meta database.
    """

    # TODO(ark): move the DB file from Docker image to S3.
    # Implement syncing procedure.
    # Update URL once the original file is available again.
    DB_PATH = "data/microarray/GEOmetadb.sqlite"
    SOURCE = "geo-meta-db"
    SOURCE_NAME = "microarray-geo"
    TECHNOLOGY = "microarray"

    def __str__(self):
        return "MicroArray GEO accession agent"

    def build_query(self) -> str:
        """Returns a query for getting GEO accessions from the local SQLite meta DB."""
        tables = [
            "SELECT *",
            "FROM gse_gpl",
            "JOIN gpl ON gse_gpl.gpl=gpl.gpl",
            "JOIN gse ON gse.gse=gse_gpl.gse",
            "GROUP BY gse_gpl.gse",
        ]

        conditions = [
            f"HAVING gse.submission_date >= '{self.since}'",
            f"gse.submission_date <= '{self.until}'",
        ]

        if self.ids:
            gpl_ids = (f"'{gpl_id}'" for gpl_id in self.ids)
            conditions.append(f"gse_gpl.gpl IN ({', '.join(gpl_ids)})")
        elif self.organism:
            conditions.append(f"lower(organism)='{self.organism.lower()}'")

        return f"{' '.join(tables)} {' AND '.join(conditions)}"

    def collect_data(self) -> Set[str]:
        """Gets new accessions from GEO database."""
        accessions = set()

        if self.ids:
            message = (
                "Getting MicroArray GEO entries by GEO platform ID(s): "
                f"{', '.join(self.ids)} for [{self.since} - {self.until}] range."
            )
        elif self.keyword:
            message = (
                f'Getting MicroArray GEO entries by "{self.keyword}" keyword '
                f"for [{self.since} - {self.until}] range."
            )
        elif self.organism:
            message = (
                f'Getting MicroArray GEO entries by "{self.organism}" organism '
                f"for [{self.since} - {self.until}] range."
            )
        else:
            return accessions

        logger.debug(message)
        accessions.update(self.fetch_data())

        return accessions

    def fetch_data(self) -> Set[str]:
        """Retrieves accessions from the GEO meta DB."""

        def match_keyword(row):
            """
            Returns True if `row` matches `self.keyword` based regex.
            Otherwise returns False.
            """
            return re_keyword.match(" ".join((str(c) for c in row if c)))

        accessions = set()

        if not os.path.exists(self.DB_PATH):
            logger.error("GEO meta database doesn't exist.")
            return accessions

        connection = sqlite3.connect(self.DB_PATH)
        connection.row_factory = sqlite3.Row
        connection.text_factory = lambda b: b.decode(errors="ignore")
        entries = connection.execute(self.build_query()).fetchall()
        connection.close()

        if self.keyword:
            re_keyword = re.compile(f".*{self.keyword}.*", re.IGNORECASE)  # Keyword regex.
            entries = filter(match_keyword, entries)

        entries = ({key.lower(): entry[key] for key in entry.keys()} for entry in entries)
        entries = set(
            (
                GatheredAccession.create_from_external_entry(entry, self.SOURCE, self.TECHNOLOGY)
                for entry in entries
            )
        )

        if self.previous_accessions:
            entries = (
                entry for entry in entries if entry.accession_code not in self.previous_accessions
            )
        accessions.update(entries)

        return accessions

    def get_ids(self) -> List[str]:
        """Returns a combined list of passed GEO platform IDs."""
        ids = set()

        if self.options["gpl_id"]:
            ids.update(self.options["gpl_id"])

        if self.options["gpl_ids_file"]:
            with open(self.options["gpl_ids_file"]) as gpl_ids_file:
                ids.update((gpl_id.strip() for gpl_id in gpl_ids_file.readlines()))

        return sorted(ids)

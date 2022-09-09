"""MicroArray (ArrayExpress, GEO) and RNA-Seq accession gathering automation.
Data sources:
  - https://www.ebi.ac.uk/biostudies/help (MicroArray ArrayExpress).
  - local SQLite meta DB from https://www.bioconductor.org/packages/release/bioc/html/GEOmetadb.html
    (MicroArray GEO).
  - https://www.ebi.ac.uk/ena/portal/api/ (RNA-Seq).
"""

import argparse
import logging
import os
import re
import sqlite3
from datetime import datetime
from http.client import RemoteDisconnected
from json.decoder import JSONDecodeError
from typing import List, Set
from urllib.parse import quote

from django.core.management.base import BaseCommand

import requests
from requests.exceptions import ConnectionError, ConnectTimeout
from retrying import retry
from urllib3.exceptions import ProtocolError

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.accession import AccessionBacklogEntry
from data_refinery_common.models.experiment import Experiment

log = get_and_configure_logger(__name__)


class Command(BaseCommand):
    """Creates agents and runs actual accession gathering."""

    RE_ACCESSION = re.compile(r"(\D+)(\d+)")
    RE_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")

    # TODO(ark): remove after upgrade to python3.8 where parser argument
    # "extend" action is directly available.
    # https://docs.python.org/3.8/library/argparse.html#action
    class ExtendAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            items = getattr(namespace, self.dest) or []
            items.extend(values)
            setattr(namespace, self.dest, items)

    def add_arguments(self, parser) -> None:
        parser.register("action", "extend", Command.ExtendAction)

        parser.add_argument(
            "--ae-id",
            action="extend",
            nargs="+",
            type=str,
            help="ArrayExpress ID(s) to use for filtering.",
        )
        parser.add_argument(
            "--ae-ids-file",
            type=str,
            help="Path to a file containing ArrayExpress ID(s) to use for filtering.",
        )
        parser.add_argument("-c", "--count", type=int, help="Number of accessions to collect.")
        parser.add_argument(
            "-d",
            "--dry-run",
            action="store_true",
            default=False,
            help="Do not write the result to the database.",
        )
        parser.add_argument(
            "-e",
            "--exclude-previous",
            action="store_true",
            default=True,
            help="Exclude previously gathered or surveyed accessions.",
        )
        parser.add_argument(
            "-ne",
            "--no-exclude-previous",
            action="store_false",
            default=False,
            dest="exclude_previous",
            help="Do not exclude previously gathered or surveyed accessions.",
        )
        parser.add_argument(
            "--gpl-id",
            action="extend",
            nargs="+",
            type=str,
            help="GEO platform ID(s) to use for filtering.",
        )
        parser.add_argument(
            "--gpl-ids-file",
            type=str,
            help="Path to a file containing GEO platform ID(s) to use for filtering.",
        )
        parser.add_argument(
            "-k",
            "--keyword",
            type=str,
            help="Keyword to use for filtering.",
        )
        parser.add_argument(
            "-m",
            "--microarray",
            action="store_true",
            default=False,
            help="Collect MicroArray accessions.",
        )
        parser.add_argument(
            "-o", "--organism", type=str, help="Organism name to use for filtering."
        )
        parser.add_argument(
            "-r",
            "--rna-seq",
            action="store_true",
            default=False,
            help="Collect RNA-Seq accessions.",
        )
        parser.add_argument(
            "-s",
            "--since",
            type=str,
            required=True,
            help="Collect accessions made public on or after this date.",
        )
        parser.add_argument(
            "--taxon-id",
            action="extend",
            nargs="+",
            type=int,
            help="Taxon ID(s) to use for filtering.",
        )
        parser.add_argument(
            "--taxon-ids-file",
            type=str,
            help="Path to a file containing taxon ID(s) to use for filtering.",
        )
        parser.add_argument(
            "-u",
            "--until",
            type=str,
            help="Collect accessions made public before or on this date.",
        )
        parser.add_argument(
            "-lv",
            "--log-verbose",
            action="store_true",
            default=False,
            help="Enable verbose log output.",
        )

    def set_verbosity_level(self, options) -> None:
        """Configures log verbosity level."""
        if options["log_verbose"]:
            log.addHandler(logging.StreamHandler())
            log.setLevel(logging.DEBUG)
        else:
            log.setLevel(logging.ERROR)

    def validate_args(self, options) -> None:
        """Validates arguments."""
        if not options["microarray"] and not options["rna_seq"]:
            exit("Either --microarray or --rna-seq must be specified.")

        errors = list()
        since = options["since"]
        until = options["until"]
        if not self.RE_DATE.match(since):
            errors.append('The -s, --since value must match "YYYY-MM-DD" format.')
        if until and not self.RE_DATE.match(until):
            errors.append('The -u, --until value must match "YYYY-MM-DD" format.')
        if since and until and since > until:
            errors.append("The -s, --since date must be earlier than -u, --until date.")

        keyword = options["keyword"]
        organism = options["organism"]
        if options["microarray"]:
            ae_id = options["ae_id"] or options["ae_ids_file"]
            gpl_id = options["gpl_id"] or options["gpl_ids_file"]
            ids = ae_id or gpl_id
            invalid_options_message = (
                "Exactly one of the keyword [-k, --keyword], organism [-o, --organism] or "
                "ArrayExpress ID(s) [--ae-id, --ae-ids-file] / GEO platform ID(s) "
                "[--gpl-id, --gpl-ids-file] must be specified."
            )
        elif options["rna_seq"]:
            taxon_id = options["taxon_id"] or options["taxon_ids_file"]
            ids = taxon_id
            invalid_options_message = (
                "Exactly one of the keyword [-k, --keyword], organism [-o, --organism] "
                "or taxon ID(s) [--taxon-id, --taxon-ids-file] must be specified."
            )

        if len([option for option in (ids, keyword, organism) if option]) != 1:
            errors.append(invalid_options_message)

        if errors:
            exit("\n".join(errors))

    def handle(self, *args, **options):
        """Runs the accession gathering process."""
        self.validate_args(options)
        self.set_verbosity_level(options)

        agents = list()
        if options["rna_seq"]:
            agents.append(RNASeqAccessionAgent(options))
        elif options["microarray"]:
            if (
                options["ae_id"]
                or options["ae_ids_file"]
                or options["keyword"]
                or options["organism"]
            ):
                agents.append(MicroArrayExpressAccessionAgent(options))
            if (
                options["gpl_id"]
                or options["gpl_ids_file"]
                or options["keyword"]
                or options["organism"]
            ):
                agents.append(MicroArrayGEOAccessionAgent(options))

        entries = set()
        for agent in agents:
            entries.update(agent.collect_data())

        entries = sorted(  # Sort the resulting list.
            (entry for entry in entries if self.RE_ACCESSION.match(entry.code)),
            key=lambda entry: (
                self.RE_ACCESSION.match(entry.code).group(1),
                int(self.RE_ACCESSION.match(entry.code).group(2)),
            ),
        )
        # Limit the number of output entries.
        entries = entries[: options["count"]] if options["count"] else entries

        if options["dry_run"]:
            if entries:
                output = "\n".join((str(entry) for entry in entries))
            else:
                output = "No accessions found."
            print(output)
        else:
            AccessionBacklogEntry.objects.bulk_create(entries)


class AccessionAgentBase:
    "Accession agent base class."

    previous_accessions = set()
    retry_params = {
        "retry_on_exception": lambda e: isinstance(
            e, (ConnectionError, ConnectTimeout, ProtocolError, RemoteDisconnected)
        ),
        "stop_max_attempt_number": 5,
        "wait_exponential_multiplier": 1000,  # Seconds.
        "wait_exponential_max": 16000,  # Seconds.
    }

    def __init__(self, options) -> None:
        """Populates args and values for major variables."""
        self.options = options
        self.count = options["count"]
        self.keyword = options["keyword"]
        self.organism = options["organism"]
        self.since = options["since"]
        self.until = options["until"] or datetime.now().strftime("%Y-%m-%d")

        self.populate_previous_accessions()

    def build_query(self):
        """Returns query/query dict depending on the accession data source."""
        raise NotImplementedError

    def collect_data(self):
        """Generates resulting entry collection."""
        raise NotImplementedError

    def fetch_data(self):
        """Fetches data from an external or local data source."""
        raise NotImplementedError

    def get_ids(self):
        """Gets IDs for query filtering depending on the accession technology."""
        raise NotImplementedError

    def populate_previous_accessions(self) -> None:
        """Populates previous accession set from a provided excluded ids file."""
        if not self.options["exclude_previous"] or self.previous_accessions:
            return

        # Gathered accessions.
        self.previous_accessions.update(
            (entry["code"] for entry in AccessionBacklogEntry.objects.values("code"))
        )

        # Surveyed accessions.
        experiments = Experiment.objects.values("accession_code", "alternate_accession_code")
        self.previous_accessions.update(
            (experiment["accession_code"] for experiment in experiments)
        )
        self.previous_accessions.update(
            (experiment["alternate_accession_code"] for experiment in experiments)
        )


class MicroArrayExpressAccessionAgent(AccessionAgentBase):
    """
    MicroArray ArrayExpress accession gathering agent. The data is fetched from
    the BioStudies database. See https://www.ebi.ac.uk/biostudies/help and
    https://www.ebi.ac.uk/biostudies/arrayexpress/help#programmatic for more
    information about the API endpoints.
    """

    def __init__(self, options) -> None:
        super().__init__(options)

        self.data_chunk_size = 100
        self.data_url = "https://www.ebi.ac.uk/biostudies/api/v1/search"
        self.ids = self.get_ids()

    def build_query(self) -> dict:
        """Returns a query dict for getting array/organism specific accessions."""
        query_dict = {
            "directsub": "true",
            "page": 1,
            "pageSize": self.data_chunk_size,
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

        log.debug(message)
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
            range_end = (params["page"] - 1) * params["pageSize"] + self.data_chunk_size
            log.debug(f"Processing entries {range_start} - {range_end}")

            response = get_response(self.data_url, params=params)
            entries = response.json().get("hits")
            if entries:
                entries = (
                    AccessionBacklogEntry.create_from_ma_ae_entry(entry) for entry in entries
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


class MicroArrayGEOAccessionAgent(AccessionAgentBase):
    """
    MicroArray GEO accession gathering agent. The data is fetched from a local
    SQLite GEO meta database.
    """

    def __init__(self, options) -> None:
        super().__init__(options)

        self.db_path = "data/microarray/GEOmetadb.sqlite"
        self.ids = self.get_ids()

    def build_query(self) -> str:
        """Returns a query for getting GEO accessions from the local SQLite meta DB."""
        tables = [
            f"SELECT *",
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

        log.debug(message)
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

        if not os.path.exists(self.db_path):
            log.error("GEO meta database doesn't exist.")
            return accessions

        connection = sqlite3.connect(self.db_path)
        connection.row_factory = sqlite3.Row
        connection.text_factory = lambda b: b.decode(errors="ignore")
        entries = connection.execute(self.build_query()).fetchall()
        connection.close()

        if self.keyword:
            re_keyword = re.compile(f".*{self.keyword}.*", re.IGNORECASE)  # Keyword regex.
            entries = filter(match_keyword, entries)

        entries = ({key.lower(): entry[key] for key in entry.keys()} for entry in entries)
        entries = set((AccessionBacklogEntry.create_from_ma_geo_entry(entry) for entry in entries))

        if self.previous_accessions:
            entries = (entry for entry in entries if entry.code not in self.previous_accessions)
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


class RNASeqAccessionAgent(AccessionAgentBase):
    """
    RNA-Seq accession gathering agent. The data is fetched from
    The European Nucleotide Archive (ENA) Portal.
    See https://www.ebi.ac.uk/ena/portal/api/ for more information about the API
    endpoints.
    """

    def __init__(self, options) -> None:
        super().__init__(options)

        self.data_chunk_size = 10000
        self.data_url = "https://www.ebi.ac.uk/ena/portal/api/search"
        self.ids = self.get_ids()

    def build_query(self, taxon_id: str = None) -> str:
        """
        Returns a query to use for getting specific taxon ID accessions.
        Some special characters must remain unquoted.
        """

        AND = " AND "
        OR = " OR "
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
            log.debug(
                f"Getting RNA-Seq entries by taxon ID(s): "
                f"{', '.join((str(idx) for idx in self.ids))} for [{self.since} - {self.until}] range."
            )
            total = len(self.ids)
            for idx, taxon_id in enumerate(self.ids):
                if self.count and len(accessions) >= self.count:
                    break

                if total > 1:
                    log.debug(f"Getting entries for taxon ID {taxon_id}, {idx + 1} of {total}.")
                accessions.update(self.fetch_data(taxon_id=taxon_id))
        elif self.keyword:
            log.debug(
                f'Getting RNA-Seq entries by "{self.keyword}" keyword '
                f"for [{self.since} - {self.until}] range."
            )
            accessions.update(self.fetch_data())
        elif self.organism:
            log.debug(
                f'Getting entries by "{self.organism}" organism '
                f"for [{self.since} - {self.until}] range."
            )
            accessions.update(self.fetch_data())

        return accessions

    def fetch_data(self, taxon_id=None) -> Set[str]:
        """
        Retrieves accessions from API search endpoint.
        The API allows to set limit to 0 (get all in one request) but we do
        it in a paginated fashion with `self.data_chunk_size` as a page size.
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
            "limit": self.data_chunk_size,
            "offset": 0,
            "query": self.build_query(taxon_id=taxon_id),
            "result": "read_study",
            "sortFields": fields,
        }

        is_done = False
        while not is_done:
            log.debug(
                f"Processing entries {data['offset'] + 1} - {data['offset'] + self.data_chunk_size}"
            )
            entries = ()
            try:
                response = get_response(self.data_url, data=data)
                entries = response.json()
                # TODO(ark): add `organism` when -o, --organism flag is used.
                entries = (
                    AccessionBacklogEntry.create_from_rnaseq_entry(entry) for entry in entries
                )
            except JSONDecodeError:
                is_done = True
            except TypeError:
                log.error(f"Couldn't get data from {self.data_url}. Response: {entries}")
            data["offset"] += self.data_chunk_size

            if self.previous_accessions:
                entries = (entry for entry in entries if entry.code not in self.previous_accessions)
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

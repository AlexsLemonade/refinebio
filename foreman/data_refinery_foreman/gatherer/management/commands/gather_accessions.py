"""MicroArray (ArrayExpress, GEO) and RNA-Seq accession gathering automation.
Data sources:
  - https://www.ebi.ac.uk/biostudies/help (MicroArray ArrayExpress).
  - local SQLite meta DB from https://www.bioconductor.org/packages/release/bioc/html/GEOmetadb.html
    (MicroArray GEO).
  - https://www.ebi.ac.uk/ena/portal/api/ (RNA-Seq).
"""

import argparse
import re

from django.core.management.base import BaseCommand
from django.db.utils import IntegrityError
from django.template.defaultfilters import pluralize

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.gathered_accession import GatheredAccession
from data_refinery_foreman.gatherer.agents.ae_agent import AEAgent
from data_refinery_foreman.gatherer.agents.geo_agent import GEOAgent
from data_refinery_foreman.gatherer.agents.rnaseq_agent import RNASeqAgent

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    """Creates agents and runs actual accession gathering."""

    DATA_AGENTS = (AEAgent, GEOAgent, RNASeqAgent)
    DATA_SOURCE_NAMES = [agent.SOURCE_NAME for agent in DATA_AGENTS]

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
            "-ne",
            "--no-exclude-previous",
            action="store_false",
            default=False,
            dest="exclude_previous",
            help="Do not exclude previously gathered or surveyed accessions.",
        )
        parser.add_argument(
            "-o", "--organism", type=str, help="Organism name to use for filtering."
        )
        parser.add_argument(
            "-s",
            "--since",
            type=str,
            required=True,
            help="Collect accessions made public on or after this date.",
        )
        parser.add_argument(
            "-src",
            "--source",
            type=str,
            action="extend",
            nargs="+",
            help="Gather accessions from selected sources.",
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

    def validate_args(self, options) -> None:
        """Validates arguments."""
        errors = list()

        RE_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
        since = options["since"]
        until = options["until"]
        if not RE_DATE.match(since):
            errors.append('The -s, --since value must match "YYYY-MM-DD" format.')
        if until and not RE_DATE.match(until):
            errors.append('The -u, --until value must match "YYYY-MM-DD" format.')
        if since and until and since > until:
            errors.append("The -s, --since date must be earlier than -u, --until date.")

        keyword = options["keyword"]
        organism = options["organism"]
        source_names = options["source"] or self.DATA_SOURCE_NAMES

        for source_name in source_names:
            if source_name in self.DATA_SOURCE_NAMES:
                continue
            errors.append(
                f"Unknown source: {source_name}. Supported sources: {', '.join(self.DATA_SOURCE_NAMES)}"
            )

        if AEAgent.SOURCE_NAME in source_names:
            ids = options["ae_id"] or options["ae_ids_file"]
            if not (ids or keyword or organism):
                errors.append(
                    (
                        "Exactly one of the keyword [-k, --keyword], organism [-o, --organism] or "
                        "ArrayExpress ID(s) [--ae-id, --ae-ids-file] must be specified for "
                        f"'{AEAgent.SOURCE_NAME}' source."
                    )
                )
        if GEOAgent.SOURCE_NAME in source_names:
            ids = options["gpl_id"] or options["gpl_ids_file"]
            if not (ids or keyword or organism):
                errors.append(
                    (
                        "Exactly one of the keyword [-k, --keyword], organism [-o, --organism] or "
                        "GEO platform ID(s) [--gpl-id, --gpl-ids-file] must be specified for "
                        f"'{GEOAgent.SOURCE_NAME}' source."
                    )
                )
        if RNASeqAgent.SOURCE_NAME in source_names:
            ids = options["taxon_id"] or options["taxon_ids_file"]
            if not (ids or keyword or organism):
                errors.append(
                    (
                        "Exactly one of the keyword [-k, --keyword], organism [-o, --organism] "
                        "or taxon ID(s) [--taxon-id, --taxon-ids-file] must be specified for "
                        f"'{RNASeqAgent.SOURCE_NAME}' source."
                    )
                )

        if errors:
            exit("\n".join(errors))

    def handle(self, *args, **options):
        """Creates agents and runs the accession gathering process."""
        self.validate_args(options)

        agents = list()
        sources_names = options["source"] or self.DATA_SOURCE_NAMES
        for cls in self.DATA_AGENTS:
            if cls.SOURCE_NAME not in sources_names:
                continue
            agents.append(cls(options))

        gathered_accessions = set()
        for agent in agents:
            agent_accessions = agent.collect_data()
            agent_accessions_count = len(agent_accessions)

            gathered_accessions.update(agent_accessions)

            logger.info(
                f"Since {options['since']} {agent} has gathered {agent_accessions_count} "
                f"new accession{pluralize(agent_accessions_count)}."
            )

        RE_ACCESSION = re.compile(r"(\D+)(\d+)")
        gathered_accessions = sorted(  # Sort the resulting list.
            (ga for ga in gathered_accessions if RE_ACCESSION.match(ga.accession_code)),
            key=lambda ga: (
                RE_ACCESSION.match(ga.accession_code).group(1),
                int(RE_ACCESSION.match(ga.accession_code).group(2)),
            ),
        )
        # Limit the number of output entries.
        gathered_accessions = (
            gathered_accessions[: options["count"]] if options["count"] else gathered_accessions
        )

        agents_count = len(agents)
        gathered_accessions_count = len(gathered_accessions)
        logger.info(
            f"Since {options['since']} {gathered_accessions_count} "
            f"new accession{pluralize(gathered_accessions_count)} "
            f"h{pluralize(gathered_accessions_count, 'as,ave')} been gathered by {agents_count} "
            f"accession agent{pluralize(agents_count)}."
        )

        if options["dry_run"]:
            if gathered_accessions:
                output = "\n".join((str(ga) for ga in gathered_accessions))
            else:
                output = "No new accessions gathered."
            exit(output)

        BATCH_SIZE = 100
        for batch_start in range(0, len(gathered_accessions), BATCH_SIZE):
            batch = gathered_accessions[batch_start : batch_start + BATCH_SIZE]
            try:
                GatheredAccession.objects.bulk_create(batch)
                added_accessions_count = len(batch)
                logger.info(
                    f"{added_accessions_count} new accession{pluralize(added_accessions_count)} "
                    f"h{pluralize(added_accessions_count, 'as,ave')} been added to the database. "
                )
                logger.info("Added accessions: '%s'" % batch)
            except IntegrityError as e:
                logger.exception(f"Could not save new accessions to the database: {e}")
                logger.info("Failed accessions: '%s'" % batch)

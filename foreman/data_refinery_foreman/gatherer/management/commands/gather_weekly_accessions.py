"""A wrapper for `gather_accessions` command which can be executed as a cron job."""

from datetime import datetime

from django.core.management import call_command
from django.core.management.base import BaseCommand

from dateutil.relativedelta import SU, relativedelta


class Command(BaseCommand):
    """Gathers accessions published since previous Sunday."""

    def handle(self, *args, **options):
        """Runs accession gathering."""
        last_sunday = (datetime.now() - relativedelta(weekday=SU(-1))).strftime("%Y-%m-%d")

        # TODO(ark): make the command process all sources within a single run.
        call_command(
            "gather_accessions",
            "--microarray",
            since=last_sunday,
            organism="homo sapiens",
        )
        call_command(
            "gather_accessions",
            "--rna-seq",
            since=last_sunday,
            taxon_ids_file="config/rna_seq_taxon_ids.txt",
        )

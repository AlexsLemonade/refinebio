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

        call_command(
            "gather_accessions",
            organism="homo sapiens",
            since=last_sunday,
            taxon_ids_file="config/rna_seq_taxon_ids.txt",
        )

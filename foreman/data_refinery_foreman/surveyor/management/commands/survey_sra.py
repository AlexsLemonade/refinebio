"""
This command will create and run survey jobs for each SRA run accession
in the range from start_accession to end_accession.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "start_accession",
            help=("An SRA run accession. Runs going from the start_accession to the "
                  "end_accession will be surveyed. If surveying a single run accession "
                  "is desired then the same accession should be passed to both start_accession "
                  "and end_accession."))
        parser.add_argument(
            "end_accession",
            help=("An SRA run accession. Runs going from the start_accession to the "
                  "end_accession will be surveyed. If surveying a single run accession "
                  "is desired then the same accession should be passed to both start_accession "
                  "and end_accession."))

    def handle(self, *args, **options):
        if options["start_accession"] is None or options["end_accession"] is None:
            logger.error("You must specify start_accession and end_accession.")
            return 1
        else:
            surveyor.survey_sra_experiments(options["start_accession"], options["end_accession"])
            return 0

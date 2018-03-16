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
            "--accession",
            type=str,
            help=("An SRA run accession. "))
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes.")
        )

    def handle(self, *args, **options):
        if options["accession"] is None and options["file"] is None:
            logger.error("You must specify accession or input file.")
            return 1
        if options["file"]:
            with open(options["file"]) as file:
                for acession in file:
                    try:
                        surveyor.survey_sra_experiment(accession.strip())
                    except Exception as e:
                        print(e)
        else:
            surveyor.survey_sra_experiment(options["accession"])
            return 0

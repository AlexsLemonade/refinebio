"""
This command will create and run survey jobs for each experiment in the
experiment_list. experiment list should be a file containing one
experiment accession code per line.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession",
            help=("An experiment accession code to survey, download, and process."))
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes.")
        )

    def handle(self, *args, **options):
        if options["accession"] is None and options['file'] is None:
            logger.error("You must specify an experiment accession or file.")
            return 1
        if options["file"]:
            with open(options["file"]) as file:
                for accession in file:
                    try:
                        surveyor.survey_ae_experiment(self.get_ae_accession(accession))
                    except Exception as e:
                        print(e)        
        else:
            surveyor.survey_ae_experiment(self.get_ae_accession(options['accession']))
            return 0

    def get_ae_accession(self, accession):
        """This allows us to support ascession codes in both
        # ArrayExpress and imported-from-GEO format."""
        if "GSE" in accession:
            accession = "E-GEOD-" + accession.split('GSE')[1] 
        return accession.strip()
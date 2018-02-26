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
            "experiment_accession",
            help=("An experiment accession code to survey, download, and process."))

    def handle(self, *args, **options):
        if options["experiment_accession"] is None:
            logger.error("You must specify an experiment accession.")
            return 1
        else:

            # This allos us to support ascession codes in both
            # ArrayExpress and imported-from-GEO format.
            accession = options["experiment_accession"]
            if "GSE" in accession:
                accession = "E-GEOD-" + accession.split('GSE')[1] 

            surveyor.survey_ae_experiment(accession)
            return 0

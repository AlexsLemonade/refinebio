"""
This command will create and run survey jobs for each experiment in the
experiment_list. experiment list should be a file containing one
experiment accession code per line.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "experiment_list",
            help=("A file containing a list of experiment accession codes to "
                  "survey, download, and process. These should be listed one "
                  "per line. Should be a path relative to the foreman "
                  "directory."))

    def handle(self, *args, **options):
        if options["experiment_list"] is None:
            logger.error("You must specify an experiment list.")
            return 0
        else:
            surveyor.survey_experiments(options["experiment_list"])

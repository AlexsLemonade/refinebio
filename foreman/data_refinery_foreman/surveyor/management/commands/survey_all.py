"""This command will create and run survey jobs for an experiment
specified by an accession. The type of survey job to run will be
determined by the pattern of the accession.
"""

import boto3
import botocore
import uuid

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)

def run_surveyor_for_accession(accession: str) -> None:
    """Chooses the correct surveyor based on the pattern of the accession"""
    if 'GSE' in accession[:3]:
        surveyor.survey_experiment(accession, "GEO")
    elif 'E-' in accession[:2]:
        surveyor.survey_experiment(accession, "ARRAY_EXPRESS")
    elif " " in accession:
        surveyor.survey_transcriptome_index(accession)
    else:
        surveyor.survey_experiment(accession, "SRA")

class Command(BaseCommand):
    def add_arguments(self, parser):

        parser.add_argument(
            "--accession",
            type=str,
            help=("An accession code to survey.")
        )

    def handle(self, *args, **options):
        if options['accession'] is None:
            logger.error("You must specify an accession.")
            return "1"

        try:
            run_surveyor_for_accession(options["accession"])
        except Exception as e:
            logger.exception(e)

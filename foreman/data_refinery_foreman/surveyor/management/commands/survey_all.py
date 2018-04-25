"""
This command will create and run survey jobs for each experiment 
in a file containing one experiment accession code per line.
"""

import boto3
import botocore
import uuid

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes. s3:// URLs are also accepted.")
        )
        parser.add_argument(
            "--offset",
            type=int,
            help=("Skip a number of lines at the beginning"),
            default=0
        )

    def handle(self, *args, **options):
        if options['file'] is None:
            logger.error("You must specify a file.")
            return 1

        if options["file"]:

            if 's3://' in options["file"]:
                bucket, key = parse_s3_url(options["file"])
                s3 = boto3.resource('s3')
                try:
                    filepath = "/tmp/input_" + str(uuid.uuid4()) + ".txt"
                    s3.Bucket(bucket).download_file(key, "/tmp/input_" + str(uuid.uuid4()) + ".txt")
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        logger.error("The remote file does not exist.")
                        raise
                    else:
                        raise
            else:
                filepath = options["file"]

            with open(filepath) as accession_file:
                for i, accession in enumerate(accession_file):
                    if i < options["offset"]:
                        continue
                    accession = accession.strip()
                    try:
                        if 'GSE' in accession[:3]:
                            surveyor.survey_geo_experiment(accession)
                        elif 'E-' in accession[:2]:
                            surveyor.survey_ae_experiment(accession)
                        else:
                            surveyor.survey_sra_experiment(accession)
                    except Exception as e:
                        logger.error(e)

"""
This command will create and run survey jobs for each SRA run accession
in the range from start_accession to end_accession.
"""

import boto3
import botocore
import uuid
import sys

from django.core.management.base import BaseCommand
from data_refinery_foreman.surveyor import surveyor
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import parse_s3_url

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
            help=("An optional file listing accession codes. s3:// URLs are also accepted.")
        )

    def handle(self, *args, **options):
        if options["accession"] is None and options["file"] is None:
            logger.error("You must specify accession or input file.")
            sys.exit(1)
        if options["file"]:
            if 's3://' in options["file"]:
                bucket, key = parse_s3_url(options["file"])
                s3 = boto3.resource('s3')
                try:
                    filepath = "/tmp/input_" + str(uuid.uuid4()) + ".txt"
                    s3.Bucket(bucket).download_file(key, filepath)
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        logger.error("The remote file does not exist.")
                    raise
            else:
                filepath = options["file"]

            with open(filepath) as file:
                for accession in file:
                    try:
                        surveyor.survey_sra_experiment(accession.strip())
                    except Exception as e:
                        print(e)
        else:
            surveyor.survey_sra_experiment(options["accession"])
            sys.exit(0)

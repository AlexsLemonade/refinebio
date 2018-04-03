"""
This command will create and run survey jobs for each experiment in the
experiment_list. experiment list should be a file containing one
experiment accession code per line.
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
            "--accession",
            help=("An experiment accession code to survey, download, and process."))
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes. s3:// URLs are also accepted.")
        )

    def handle(self, *args, **options):
        if options["accession"] is None and options['file'] is None:
            logger.error("You must specify an experiment accession or file.")
            return 1
        if options["file"]:

            if 's3://' in options["file"]:
                bucket, key = parse_s3_url(file)
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

            with open(filepath) as file:
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
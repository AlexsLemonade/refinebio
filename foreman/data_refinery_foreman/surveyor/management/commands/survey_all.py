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
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)


def run_surveyor_for_accession(accession: str) -> None:
    """Chooses the correct surveyor based on the pattern of the accession"""
    if "GSE" in accession[:3]:
        surveyor.survey_experiment(accession, "GEO")
    elif "E-" in accession[:2]:
        surveyor.survey_experiment(accession, "ARRAY_EXPRESS")
    elif " " in accession:
        args = accession.split(",")
        # Allow organism to be unspecified so we survey the entire division.
        organism = args[0] if len(args[0]) > 0 else None
        if len(args) > 1:
            division = args[1].strip()
        else:
            division = "Ensembl"
        surveyor.survey_transcriptome_index(organism, division)
    else:
        surveyor.survey_experiment(accession, "SRA")


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help=("""An optional file listing accession codes. s3:// URLs are also accepted."""),
        )
        parser.add_argument("--accession", type=str, help=("An accession code to survey."))
        parser.add_argument(
            "--offset", type=int, help=("Skip a number of lines at the beginning"), default=0,
        )
        parser.add_argument(
            "--job-id", type=int, help=("An ID of a SurveyJob to execute"), default=None
        )

    def handle(self, *args, **options):
        if options["file"] is None and options["accession"] is None and options["job_id"] is None:
            logger.error("You must specify an accession or file or job ID.")
            return "1"

        if options["file"]:
            if "s3://" in options["file"]:
                bucket, key = parse_s3_url(options["file"])
                s3 = boto3.resource("s3")
                try:
                    filepath = "/tmp/input_" + str(uuid.uuid4()) + ".txt"
                    s3.Bucket(bucket).download_file(key, filepath)
                except botocore.exceptions.ClientError as e:
                    if e.response["Error"]["Code"] == "404":
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
                        run_surveyor_for_accession(accession)
                    except Exception as e:
                        logger.exception(e)

        if options["accession"]:
            accession = options["accession"]
            try:
                run_surveyor_for_accession(accession)
            except Exception as e:
                logger.exception(e)

        if options["job_id"]:
            job_id = options["job_id"]
            try:
                survey_job = SurveyJob.objects.get(id=job_id)
                surveyor.run_job(survey_job)
            except Exception as e:
                logger.exception(e)

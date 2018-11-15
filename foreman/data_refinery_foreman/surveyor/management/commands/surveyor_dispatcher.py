"""This command will create and run survey jobs for each experiment
in the experiment_list. experiment list should be a file containing
one experiment accession code per line.
"""

import boto3
import botocore
import nomad
import uuid

from django.core.management.base import BaseCommand
from nomad.api.exceptions import URLNotFoundNomadException

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_lookup import SurveyJobTypes
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_common.utils import parse_s3_url, get_env_variable
from data_refinery_foreman.surveyor import surveyor


logger = get_and_configure_logger(__name__)

SURVEYOR_JOB_NAME = "SURVEYOR_256"

def set_source_type_for_accession(survey_job, accession: str) -> None:
    """Type a surveyor based on accession structure"""
    if 'GSE' in accession[:3]:
        survey_job.source_type = "GEO"
        survey_job.save()
        return
    elif 'E-' in accession[:2]:
        survey_job.source_type = "ARRAY_EXPRESS"
        survey_job.save()
        return
    elif " " in accession:

        survey_job.source_type = "TRANSCRIPTOME_INDEX"
        survey_job.save()

        args = accession.split(",")
        # Allow organism to be unspecified so we survey the entire division.
        organism_name = args[0] if len(args[0]) > 0 else None
        if len(args) > 1:
            ensembl_division = args[1].strip()
        else:
            ensembl_division = "Ensembl"

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="ensembl_division",
                                           value=ensembl_division)
        key_value_pair.save()
        if organism_name:
            key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                               key="organism_name",
                                               value=organism_name)
            key_value_pair.save()

        return
    else:
        survey_job.source_type = "SRA"
        survey_job.save()
        return

def queue_surveyor_for_accession(accession: str) -> None:
    """Dispatches a surveyor job for the accession code."""
    # Start at 256MB of RAM for surveyor jobs.
    survey_job = SurveyJob(ram_amount=256)

    set_source_type_for_accession(survey_job, accession)

    key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                       key="experiment_accession_code",
                                       value=accession)
    key_value_pair.save()

    try:
        send_job(SurveyJobTypes.SURVEYOR, survey_job)
    except:
        # If we can't dispatch this, then let the foreman retry it late.
        pass

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help=("""An file listing accession codes. s3:// URLs are also accepted."

Note: One entry per line, GSE* entries survey GEO, E-GEO-* entries survey ArrayExpress.
""")
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
            return "1"

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
                    raise
        else:
            filepath = options["file"]

        with open(filepath) as accession_file:
            for i, accession in enumerate(accession_file):
                if i < options["offset"]:
                    continue
                accession = accession.strip()
                try:
                    queue_surveyor_for_accession(accession)
                except Exception as e:
                    logger.exception(e)

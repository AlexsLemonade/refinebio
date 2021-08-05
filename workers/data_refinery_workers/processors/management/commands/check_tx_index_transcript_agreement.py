import csv
import itertools
import json
import os.path
import tarfile

from django.conf import settings
from django.core.management.base import BaseCommand

import boto3
import requests
from botocore.client import Config

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Organism, OrganismIndex

logger = get_and_configure_logger(__name__)

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))
S3_BUCKET = "data-refinery-test-assets"


def get_s3_key(organism, index_type):
    return f"/tx-index-test-{settings.USER}-{settings.STAGE}-{organism.name}-{index_type}.json"


def get_old_transcripts(organism, index_type):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Tried to get old transcripts when we are not running in the cloud.")
        return None

    key = get_s3_key(organism, index_type)
    try:
        response = S3.get_object(Bucket=S3_BUCKET, Key=key)

        return set(json.loads(response["Body"].read()))

    except boto3.S3.Client.exceptions.NoSuchKey:
        return None

    except Exception:
        logger.exception(
            "Error downloading previous test run from S3", s3_key=key, s3_bucket=S3_BUCKET
        )
        raise

    return None


def upload_new_transcripts(organism, index_type, transcripts):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Tried to upload new transcripts when we are not running in the cloud.")
        return

    key = get_s3_key(organism, index_type)
    try:
        S3.put_object(Body=json.dumps(list(transcripts)), Bucket=S3_BUCKET, Key=key)

    except Exception:
        logger.exception("Error uploading new transcripts to S3", s3_key=key, s3_bucket=S3_BUCKET)
        raise
    pass


def signal_failure(bad_indices):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Cannot message Slack when we are not running in the cloud")
        return
    elif not settings.ENGAGEMENTBOT_WEBHOOK:
        logger.error("Cannot message slack because we don't have a valid webhook")
        return

    requests.post(
        settings.ENGAGEMENTBOT_WEBHOOK,
        json={
            "attachments": [
                {
                    "fallback": "Transcriptome Index Test Failed",
                    "title": "Transcriptome Index Test Failed",
                    "color": "#db3b28",
                    "text": "Some transcriptome indices did not have enough overlap with the previous run",
                    "fields": [
                        {
                            "title": "Failing Indixes",
                            "value": ", ".join(f"{o.name} {t}" for o, t in bad_indices),
                        }
                    ],
                    "footer": "Refine.bio",
                    "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                }
            ],
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )


def get_current_transcripts(organism, index_type):
    index_object = (
        OrganismIndex.objects.filter(organism=organism, index_type=index_type)
        .order_by("-created_at")
        .first()
    )
    if not index_object:
        return None

    index_tarball = index_object.get_computed_file().sync_from_s3()

    with tarfile.open(index_tarball, "r:gz") as index_archive:
        with index_archive.open("genes_to_transcripts.txt", "r") as f:
            reader = csv.reader(f, delimiter="\t")
            transcripts = {t for _, t in reader}

    os.path.remove(index_tarball)
    return transcripts


def check_tx_index_transcript_agreement():
    bad_indices = []

    for organism, index_type in itertools.product(Organism.objects.all(), ["short", "long"]):
        current_transcripts = get_current_transcripts(organism, index_type)
        if current_transcripts is None:
            continue

        old_transcripts = get_old_transcripts(organism, index_type)
        if old_transcripts is None:
            logger.info(
                f"Skipping jaccard index test for {organism.name} for {index_type} because there is no old data"
            )
        else:
            jaccard_index = len(current_transcripts & old_transcripts) / len(
                current_transcripts | old_transcripts
            )

            if jaccard_index < 0.95:
                bad_indices.append((organism, index_type))

        upload_new_transcripts(organism, index_type, current_transcripts)

    if len(bad_indices) != 0:
        signal_failure(bad_indices)


class Command(BaseCommand):
    def handle(self, *args, **options):
        try:
            check_tx_index_transcript_agreement()
        except Exception as e:
            requests.post(
                settings.ENGAGEMENTBOT_WEBHOOK,
                json={
                    "attachments": [
                        {
                            "fallback": "Exception raised during Transcriptome Index Test",
                            "title": "Exception raised during Transcriptome Index Test",
                            "color": "#db3b28",
                            "text": f"```\n{str(e)}\n```",
                            "footer": "Refine.bio",
                            "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                        }
                    ],
                },
                headers={"Content-Type": "application/json"},
                timeout=10,
            )

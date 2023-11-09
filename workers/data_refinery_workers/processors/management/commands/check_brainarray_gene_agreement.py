import csv
import io
import json
from typing import Iterator, List, Tuple

from django.conf import settings
from django.core.management.base import BaseCommand

import boto3
import botocore
import requests
import rpy2.robjects as ro
from botocore.client import Config
from botocore.exceptions import ClientError
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))
S3_BUCKET = f"data-refinery-s3-{settings.USER}-{settings.STAGE}"


def brainarray_platforms() -> Iterator[Tuple[str, str]]:
    with open("/home/user/r_ensg_probe_pkgs.txt", "r") as f:
        for row in csv.reader(f, delimiter="\t"):
            platform = row[0]

            url = row[1]
            package_filename = url.split("/")[-1]
            package = package_filename.split("_")[0]

            yield (platform, package)


def probe_set_name_starts_with_one_of(gene_prefixes):
    """Returns a string that can be used with pandas.DataFrame.query"""
    return " | ".join(
        f"`Probe.Set.Name`.str.startswith('{prefix}', na=False)" for prefix in gene_prefixes
    )


def get_genes_for_package(package: str) -> List[str]:
    ro.r(f'suppressPackageStartupMessages(library("{package}", character.only=TRUE))')
    brainarray_df = ro.r(f'as.data.frame(get("{package}"), stringsAsFactors = FALSE)')
    with localconverter(ro.default_converter + pandas2ri.converter):
        brainarray_df = ro.conversion.rpy2py(brainarray_df)

    gene_prefixes = ["ENS"]

    # Some species have their own special gene prefixes
    if package.startswith("celegans") or package.startswith("elegene"):
        gene_prefixes.append("WBGene")
    elif (
        package.startswith("drosgenome")
        or package.startswith("drosophila")
        or package.startswith("drogene")
    ):
        gene_prefixes.append("FBgn")
    elif package.startswith("yeast") or package.startswith("ygs"):
        # For some reason, yeast genes just have a Y prefix
        gene_prefixes.append("Y")

    genes = set(
        brainarray_df.query(probe_set_name_starts_with_one_of(gene_prefixes))
        .assign(gene=lambda df: df["Probe.Set.Name"].str.replace("_at", ""))
        .gene
    )

    return genes


def get_s3_key(platform):
    return f"brainarray-test/{platform}.json"


def get_old_platform_genes(platform):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Tried to get old platform genes when we are not running in the cloud.")
        return None

    key = get_s3_key(platform)
    try:
        response = S3.get_object(Bucket=S3_BUCKET, Key=key)
        return set(json.loads(response["Body"].read()))

    except S3.exceptions.NoSuchKey:
        return None

    # For some reason, with our current authentication scheme this gets thrown
    # instead of NoSuchKey when the key doesn't exist.
    except botocore.exceptions.ClientError:
        return None

    except Exception:
        logger.exception(
            "Error downloading previous test run from S3", s3_key=key, s3_bucket=S3_BUCKET
        )
        raise


def upload_new_platform_genes(platform, platform_genes):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Tried to upload new platform genes when we are not running in the cloud.")
        return

    key = get_s3_key(platform)
    try:
        data = json.dumps(list(platform_genes))
        S3.upload_fileobj(io.BytesIO(data.encode()), S3_BUCKET, key)

    except Exception:
        logger.exception(
            "Error uploading new platform genes to S3", s3_key=key, s3_bucket=S3_BUCKET
        )
        raise
    pass


def signal_failure(bad_platforms):
    if not settings.RUNNING_IN_CLOUD or not settings.SLACK_WEBHOOK_URL:
        logger.info(
            "Some Brainarray platforms did not have enough overlap with the previous run",
            bad_platforms=bad_platforms,
        )
        return

    requests.post(
        settings.SLACK_WEBHOOK_URL,
        json={
            "attachments": [
                {
                    "fallback": "Brianarray Test Failed",
                    "title": "Brainarray Test Failed",
                    "color": "#db3b28",
                    "text": "Some Brainarray platforms did not have enough overlap with the previous run",
                    "fields": [{"title": "Failing Platforms", "value": ", ".join(bad_platforms)}],
                    "footer": "Refine.bio",
                    "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                }
            ],
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )


def check_brainarray_gene_agreement():
    bad_platforms = []

    for platform, package in brainarray_platforms():
        logger.info(f"Checking platform {platform}...")

        platform_genes = get_genes_for_package(package)

        old_platform_genes = get_old_platform_genes(platform)
        if old_platform_genes is None:
            logger.info(f"Skipping jaccard index test for {platform} because there is no old data")
        else:
            jaccard_index = len(platform_genes & old_platform_genes) / len(
                platform_genes | old_platform_genes
            )

            if jaccard_index < 0.95:
                logger.error(
                    "Found a platform with bad overlap",
                    platform=platform,
                    jaccard_index=jaccard_index,
                )
                bad_platforms.append(platform)

                # Don't overwrite the old data if we have bad overlap, so we can run this again
                continue
            else:
                logger.info("Good overlap", platform=platform, jaccard_index=jaccard_index)

        upload_new_platform_genes(platform, platform_genes)

    if len(bad_platforms) != 0:
        signal_failure(bad_platforms)


class Command(BaseCommand):
    def handle(self, *args, **options):
        try:
            check_brainarray_gene_agreement()
        except Exception as e:
            if settings.SLACK_WEBHOOK_URL:
                requests.post(
                    settings.SLACK_WEBHOOK_URL,
                    json={
                        "attachments": [
                            {
                                "fallback": "Exception raised during Brainarray Test",
                                "title": "Exception raised during Brainarray Test",
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
            else:
                logger.exception("Exception was raised during Brainarray Test")

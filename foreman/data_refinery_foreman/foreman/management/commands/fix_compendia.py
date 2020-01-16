""" This management command will transfer all compendias in the bucket
`data-refinery-s3-circleci-prod` to `data-refinery-s3-compendia-circleci-prod`
It should be a one off fix for https://github.com/AlexsLemonade/refinebio/issues/2072
We considered adding it in a migration but we didn't liked the idea of making
S3 operations there.
"""

from django.conf import settings
from django.core.management.base import BaseCommand

import boto3
from botocore.client import Config

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputedFile
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)
S3_COMPENDIA_BUCKET_NAME = get_env_variable("S3_COMPENDIA_BUCKET_NAME", "data-refinery")

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))


def transfer_computed_file(computed_file):
    old_s3_bucket = computed_file.s3_bucket
    s3_key = computed_file.s3_key

    # 1. copy to new location
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3.html#S3.Client.copy
    copy_source = {"Bucket": old_s3_bucket, "Key": s3_key}
    # keep the same key, just different
    S3.copy(
        copy_source,
        S3_COMPENDIA_BUCKET_NAME,
        s3_key,
        ExtraArgs={"ACL": "public-read", "StorageClass": "STANDARD_IA"},
    )

    # 2. update computed file
    computed_file.s3_bucket = S3_COMPENDIA_BUCKET_NAME
    computed_file.save()

    # 3. Remove file from original bucket
    S3.delete_object(Bucket=old_s3_bucket, Key=s3_key)


class Command(BaseCommand):
    def handle(self, *args, **options):
        compendia_files = ComputedFile.objects.filter(
            result__compendium_result__isnull=False
        ).exclude(s3_bucket=S3_COMPENDIA_BUCKET_NAME)

        for file in compendia_files:
            transfer_computed_file(file)

            logger.info("Transferred compendia file " + file.filename)

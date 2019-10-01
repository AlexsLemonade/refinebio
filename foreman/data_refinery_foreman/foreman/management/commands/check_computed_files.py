from multiprocessing import Pool

import boto3
from botocore.client import Config
from botocore.errorfactory import ClientError
from django.core.management.base import BaseCommand
from django.db import connections

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputedFile
from data_refinery_common.utils import queryset_page_iterator

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client('s3', config=Config(signature_version='s3v4'))

logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000

def check_page(computed_file_page):
    missing_file_ids = []
    for computed_file in computed_file_page:
        # check that file is present in S3, no need to download the entire object
        # https://stackoverflow.com/a/38376288/763705
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/s3.html#S3.Client.head_object
        try:
            S3.head_object(Bucket=computed_file.s3_bucket, Key=computed_file.s3_key)
        except ClientError:
            # Not found
            logger.info('Computed file not found on S3 - will remove S3 fields.', computed_file=computed_file, filename=computed_file.filename)
            missing_file_ids.append(computed_file.pk)

    if missing_file_ids:
        # Update all computed files in one query
        logger.info('Checked one page of computed files, found %i were missing in S3 clearing s3_key and s3_bucket now.' % len(missing_file_ids))
        ComputedFile.objects.filter(id__in=missing_file_ids).update(s3_key=None, s3_bucket=None)
    else:
        logger.info('Checked one page of computed files, all of them were in S3.')

class Command(BaseCommand):

    def handle(self, *args, **options):
        """ We found some computed files that seem to be in S3 but they don't exist for some reason.
        This command will check all of them and ensure they exists in S3. For those that don't exist it
        will clear their `s3_key`/`s3_bucket` fields.
        ref https://github.com/AlexsLemonade/refinebio/issues/1696
        """
        computed_files_queryset = ComputedFile.objects.filter(s3_key__isnull=False, s3_bucket__isnull=False)
        results = queryset_page_iterator(computed_files_queryset, PAGE_SIZE)

        connections.close_all()
        with Pool(processes=8) as pool:
            pool.apply(check_page, results)

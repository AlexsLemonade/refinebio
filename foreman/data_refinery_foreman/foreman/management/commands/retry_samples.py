import random
import sys
import time
import re
from typing import Dict, List

from django.core.management.base import BaseCommand
from django.db.models import OuterRef, Subquery, Count
from dateutil.parser import parse as parse_date
from data_refinery_common.models import ProcessorJob, Sample, ComputedFile
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_management import create_downloader_job
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator


logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000


def requeue_samples(eligible_samples):
    paginator = Paginator(eligible_samples, PAGE_SIZE)
    page = paginator.page()

    creation_count = 0
    while True:
        for sample in page.object_list:
            if create_downloader_job(sample.original_files.all(), force=True):
                creation_count += 1

        if not page.has_next():
            break
        else:
            page = paginator.page(page.next_page_number())

        logger.info("Creating new downloader jobs. %d so far", creation_count)

        # 2000 samples queued up every five minutes should be fast
        # enough and also not thrash the DB.
        time.sleep(60 * 5)

    return creation_count


def retry_by_source_database(source_database):
    sra_samples = (
        Sample.objects.filter(source_database=source_database, is_processed=False)
        .annotate(num_computed_files=Count("computed_files"))
        .filter(num_computed_files__gt=0)
        .prefetch_related("original_files")
    )

    creation_count = requeue_samples(sra_samples)
    logger.info(
        "Created %d new downloader jobs because their samples lacked computed files.",
        creation_count,
    )


def retry_by_accession_codes(accession_codes_param):
    accession_codes = accession_codes_param.split(",")
    eligible_samples = Sample.objects.filter(accession_code__in=accession_codes).prefetch_related(
        "original_files"
    )
    total_samples_queued = requeue_samples(eligible_samples)
    logger.info(
        "Re-queued %d samples with accession codes %s.", total_samples_queued, accession_codes_param
    )


def retry_by_regex(pattern, last_modified=None):
    """ Finds the samples where the failure reason of the latest processor job that was applied to them
    matches the given regex
    https://docs.djangoproject.com/en/dev/ref/models/querysets/#regex
    """
    if last_modified:
        processor_jobs = ProcessorJob.objects.filter(last_modified__gt=last_modified)
    else:
        processor_jobs = ProcessorJob.objects

    latest_processor_job_for_sample = processor_jobs.filter(
        start_time__isnull=False, original_files__samples=OuterRef("id")
    ).order_by("-start_time")

    eligible_samples = (
        Sample.objects.filter(is_processed=False)
        .annotate(
            failure_reason=Subquery(latest_processor_job_for_sample.values("failure_reason")[:1])
        )
        .filter(failure_reason__regex=pattern)
        .prefetch_related("original_files")
    )

    total_samples_queued = requeue_samples(eligible_samples)

    logger.info(
        "Re-queued %d samples that had failed with the pattern %s.", total_samples_queued, pattern
    )


def retry_computed_files_not_uploaded(last_modified=None):
    if last_modified:
        computed_files = ComputedFile.objects.filter(last_modified__gt=last_modified)
    else:
        computed_files = ComputedFile.objects

    samples_with_computed_files = computed_files.filter(
        s3_bucket__isnull=True, samples__isnull=False
    ).values_list("samples", flat=True)

    samples_with_results = computed_files.filter(
        s3_bucket__isnull=True, result__samples__isnull=False
    ).values_list("result__samples", flat=True)

    sample_ids = set(samples_with_results) | set(samples_with_computed_files)

    eligible_samples = Sample.objects.filter(id__in=sample_ids)

    total_samples_queued = requeue_samples(eligible_samples)
    logger.info(
        "Re-queued %d samples that were associated with computed files that were not uploaded to S3",
        total_samples_queued,
    )


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--source-database",
            type=str,
            help=(
                "The name of a source database, such as Array Express, GEO, or SRA."
                "All samples from this source database will have downloader "
                "jobs requeued for them."
            ),
        )
        parser.add_argument(
            "--accession-codes",
            type=str,
            help=("Comma separated sample accession codes that need to be requeued."),
        )
        # https://docs.djangoproject.com/en/dev/ref/models/querysets/#regex
        parser.add_argument(
            "--failure-regex",
            type=str,
            help=(
                "Re-queues all samples where the latest processor job failure reason matches this regex."
            ),
        )
        # https://github.com/AlexsLemonade/refinebio/issues/1627
        parser.add_argument(
            "--computed-files-not-uploaded",
            action="store_true",
            help="Finds the samples that are associated with computed files that were not uploaded to S3 and requeue them",
        )

        parser.add_argument(
            "--last-modified",
            type=parse_date,
            help="Only recreated computed files or retry processor jobs modified after this date",
        )

    def handle(self, *args, **options):
        """ Re-queues all unprocessed RNA-Seq samples for an organism. """
        if options["source_database"]:
            retry_by_source_database(options["source_database"])

        if options["accession_codes"]:
            # --accession-codes="GSM12323,GSM98586"
            retry_by_accession_codes(options["accession_codes"])

        if options["failure_regex"]:
            # Examples
            # --failure-regex="ProcessorJob has already completed .*"
            # --failure-regex="Encountered error in R code while running AFFY_TO_PCL pipeline .*"
            retry_by_regex(options["failure_regex"], options.get("last_modified", None))

        if options["computed_files_not_uploaded"]:
            retry_computed_files_not_uploaded(options.get("last_modified", None))

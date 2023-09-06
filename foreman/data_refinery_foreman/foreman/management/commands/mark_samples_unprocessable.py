"""This command will search through all unprocessed samples and check
their most recent processor job's failure reason to see if it was a
type that we know indicates that the sample was either corrupted or of
low enough quality that we cannot hope to ever process it. This will
be stored on the sample's is_unable_to_be_processed field.
"""

from django.core.management.base import BaseCommand

from data_refinery_common.job_management import is_job_unprocessable
from data_refinery_common.models import Sample
from data_refinery_common.performant_pagination.pagination import PAGE_SIZE
from data_refinery_common.utils import queryset_page_iterator


def mark_unprocessable_samples():
    """This function performs the function explained at the head of this file.

    It does so by following this general strategy:
    1. Find samples that have is_processed=True
    2. Find the most recent processor job for that sample.
    3. Check that processor job's failure reason and mark
       the sample accordingly.

    Along the way it will also populate the following fields:
    * last_processor_job
    * last_downloader_job
    * most_recent_smashable_file
    * most_recent_quant_file
    """
    unprocessed_samples = Sample.objects.filter(is_processed=False).prefetch_related(
        "original_files__processor_jobs"
    )

    results = queryset_page_iterator(unprocessed_samples, PAGE_SIZE)

    for page in results:
        for sample in page:
            sample.last_processor_job = sample.get_most_recent_processor_job()
            sample.last_downloader_job = sample.get_most_recent_downloader_job()
            sample.most_recent_smashable_file = sample.get_most_recent_smashable_result_file()
            sample.most_recent_quant_file = sample.get_most_recent_quant_sf_file()

            if not sample.is_processed and sample.last_processor_job:
                sample.is_unable_to_be_processed = is_job_unprocessable(sample.last_processor_job)

            sample.save()


class Command(BaseCommand):
    def handle(self, *args, **options):
        """This is just the entrypoint for this management command.

        All of its work is done in a separate function because that
        makes it much easier to test."""
        mark_unprocessable_samples()

"""
This command finds Samples that were created and didn't spawned any downloader jobs.
We tried to debug the reasons why this happened on
https://github.com/alexslemonade/refinebio/issues/1391
without any luck.
"""

from django.core.management.base import BaseCommand
from django.db.models import Count
from dateutil.parser import parse as parse_date
import time

from data_refinery_common.models import Sample
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_management import create_downloader_job
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator

from data_refinery_common import job_lookup
from data_refinery_common.utils import get_supported_microarray_platforms, get_supported_rnaseq_platforms

logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--created-after',
                            type=parse_date,
                            help='Only recreate jobs created after this date')

    def handle(self, *args, **options):
        """ Requeues downloader jobs for samples that haven't been processed and their original files
        have no no downloader jobs associated with them
        """
        supported_microarray_platforms = [x['platform_accession'] for x in get_supported_microarray_platforms()]
        supported_rnaseq_platforms = [x.replace(' ', '') for x in get_supported_rnaseq_platforms()]
        all_supported_platforms = supported_microarray_platforms + supported_rnaseq_platforms # https://www.postgresql.org/docs/9.1/functions-array.html

        # Ensure selected samples have valid platforms
        samples_without_downloader = Sample.objects.all()\
                                                   .filter(platform_accession_code__in=all_supported_platforms)\
                                                   .annotate(original_files_count=Count('original_files'), downloader_job_count=Count('original_files__downloader_jobs'))\
                                                   .filter(is_processed=False, original_files_count__gt=0, downloader_job_count=0)\

        if options.get('created_after', None):
            samples_without_downloader = samples_without_downloader.filter(created_at__gt=options['created_after'])

        samples_without_downloader = samples_without_downloader.prefetch_related("original_files")

        logger.info("Found %d samples without downloader jobs, starting to create them now.", samples_without_downloader.count())

        paginator = Paginator(samples_without_downloader, PAGE_SIZE)
        page = paginator.page()

        while True:
            for sample in page.object_list:
                logger.debug("Creating downloader job for a sample.", sample=sample.accession_code)
                create_downloader_job(sample.original_files.all())

            logger.info("Created %d new downloader jobs because their samples didn't have any.", PAGE_SIZE)

            if not page.has_next():
                break

            page = paginator.page(page.next_page_number())

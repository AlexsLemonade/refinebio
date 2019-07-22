"""
This command finds Samples that were created and didn't spawned any downloader jobs.
We tried to debug the reasons why this happened on
https://github.com/alexslemonade/refinebio/issues/1391
without any luck. 
"""

from django.core.management.base import BaseCommand

from data_refinery_common.models import (
  Sample,
  DownloaderJob,
  ProcessorJob,
  OriginalFile,
  DownloaderJobOriginalFileAssociation,
  ProcessorJobOriginalFileAssociation
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_management import create_downloader_job
from data_refinery_foreman.foreman.performant_pagination.pagination import PerformantPaginator as Paginator

logger = get_and_configure_logger(__name__)

PAGE_SIZE=2000

class Command(BaseCommand):
    def handle(self, *args, **options):
      """ Requeues downloader jobs for samples that haven't been processed and their original files 
      have no no downloader jobs associated with them
      """
      samples_without_downloader = Sample.objects.all()\
          .annotate(original_files_count=Count('original_files'), downloader_job_count=Count('original_files__downloader_jobs'))\
          .filter(is_processed=False, original_files_count__gt=0, downloader_job_count=0)\
          .prefetch_related(
            "original_files"
          )

      logger.info("Found %d samples without downloader jobs, starting to create them now.", samples_without_downloader.count())
      
      paginator = Paginator(samples_without_downloader, PAGE_SIZE)
      page = paginator.page()
      page_count = 0

      while True:
          for sample in page.object_list:
            logger.debug("Creating downloader job for a sample.", sample=sample.accession_code)
            create_downloader_job(sample.original_files.all())

          logger.info("Created %d new downloader jobs because their samples didn't have any.", PAGE_SIZE)

          if not page.has_next():
              break
              
          page = paginator.page(page.next_page_number())

          # 2000 samples queued up every five minutes should be fast
          # enough and also not thrash the DB.
          time.sleep(60 * 5)


import random
import sys
import time
from typing import Dict, List

from django.core.management.base import BaseCommand

from data_refinery_common.models import (
    Sample,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_management import create_downloader_job
from data_refinery_foreman.foreman.performant_pagination.pagination import PerformantPaginator as Paginator


logger = get_and_configure_logger(__name__)

PAGE_SIZE=2000

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--source-database",
            type=str,
            help=("The name of a source database, such as Array Express, GEO, or SRA."
                  "All samples from this source database will have downloader "
                  "jobs requeued for them.")
        )

    def handle(self, *args, **options):
        """Requeues all unprocessed RNA-Seq samples for an organism.
        """
        if options["source_database"] is None:
            logger.error("You must specify a source-database.")
            sys.exit(1)
        else:
            source_database = options["source_database"]

        sra_samples = Sample.objects.filter(
            source_database=source_database
        ).prefetch_related(
            "computed_files",
            "original_files"
        )

        paginator = Paginator(sra_samples, PAGE_SIZE)
        page = paginator.page()
        page_count = 0

        creation_count_loop = 0
        creation_count_since_sleep = 0
        while True:
            for sample in page.object_list:
                if sample.computed_files.count() == 0:
                    logger.debug("Creating downloader job for a sample.",
                                 sample=sample.accession_code)
                    if create_downloader_job(sample.original_files.all(), force=True):
                        creation_count_loop += 1

            logger.info(
                "Created %d new downloader jobs because their samples lacked computed files.",
                creation_count_loop
            )

            if not page.has_next():
                break
            else:
                page = paginator.page(page.next_page_number())

            creation_count_since_sleep += creation_count_loop
            creation_count_loop = 0

            # 1000 samples queued up every five minutes should be fast
            # enough and also not thrash the DB.
            if creation_count_since_sleep >= 1000:
                time.sleep(60 * 5)
                creation_count_since_sleep = 0

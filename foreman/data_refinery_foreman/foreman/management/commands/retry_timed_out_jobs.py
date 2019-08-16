"""This command will slowly retry Salmon jobs that timed out.
This is now necessary because samples with unmated reads will no longer cause
us to time out. It will only queue 300 an hour so as to not overload ENA.
"""

import time
from typing import List

from django.core.management.base import BaseCommand

from data_refinery_foreman.foreman.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import  ProcessorJob


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        timed_out_jobs = ProcessorJob.objects.filter(
            success='f',
            failure_reason='Salmon timed out because it failed to complete within 3 hours.',
            retried_job_id__isnull=True # Only get jobs that weren't retried.
        )

        total = 0
        i = 0
        paginator = Paginator(timed_out_jobs, 1000)
        page = paginator.page()
        while page:
            for processor_job in page.object_list:
                # Reset the job so it'll start at 12288 RAM (since it'll go up from here).
                processor_job.ram_amount = 8192
                # And so it'll be retried twice more.
                processor_job.num_retries = 0
                processor_job.retried = False

                # We don't actually have to send this off to Nomad ourselves.
                # The Foreman will find it and requeue it for us!
                processor_job.save()

                if page.has_next():
                    page = paginator.page(page.next_page_number())
                else:
                    break

                total += 1
                i += 1
                # Only queue 300 of these an hour so we don't overload ENA.
                if i == 300:
                    logger.info("Requeued 300 more jobs (total %d). Sleeping for 1 hour.", total)
                    time.sleep(60*60)
                    i = 0

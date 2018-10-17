import sys

from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ProcessorJob,
)
from data_refinery_workers.processors.janitor import run_janitor
logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        return

    def handle(self, *args, **options):

        pj = ProcessorJob()
        pj.pipeline_applied = "JANITOR"
        pj.save()

        final_context = run_janitor(pj.pk)

        print("Removed: ")
        for item in final_context['deleted_items']:
            print('\t - ' + item)

        sys.exit(0)

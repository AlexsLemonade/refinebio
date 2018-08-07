import sys

from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job

from data_refinery_common.models import (
    OriginalFile
)

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):

    def handle(self, *args, **options):
        """ Deletes local original files
        """

    for original_file in OriginalFile.objects.all():
        try:
            original_file.delete_local_file()
        except Exception as e:
            print("[" + str(original_file) + "]: " + str(e))
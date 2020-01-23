import sys

from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OriginalFile

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        """
        Deletes local original files
        """
        if (
            not input("You're about to delete all of the local files! Proceed? (y/n): ")
            .lower()
            .strip()[0]
            == "y"
        ):
            sys.exit(1)

        for original_file in OriginalFile.objects.all():
            try:
                original_file.delete_local_file()
            except Exception as e:
                print("[" + str(original_file) + "]: " + str(e))

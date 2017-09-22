from django.core.management.base import BaseCommand
from data_refinery_workers.downloaders.array_express import download_array_express


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--job-id",
            type=int,
            help=("The downloader job's ID."))

    def handle(self, *args, **options):
        if options["job_id"] is None:
            logger.error("You must specify a job ID.")
            return 1
        else:
            download_array_express(options["job_id"])
            return 0

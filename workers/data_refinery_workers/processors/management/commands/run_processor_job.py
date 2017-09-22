from django.core.management.base import BaseCommand
from data_refinery_workers.processors.array_express import affy_to_pcl


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--job-id",
            type=int,
            help=("The processor job's ID."))

    def handle(self, *args, **options):
        if options["job_id"] is None:
            logger.error("You must specify a job ID.")
            return 1
        else:
            affy_to_pcl(options["job_id"])
            return 0

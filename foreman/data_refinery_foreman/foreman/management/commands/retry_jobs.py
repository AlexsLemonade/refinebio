"""
This command will run the Foreman's main function monitor_jobs.
This will cause the Foreman to check for a number of different
failures for both the DownloaderJobs and ProcessorJobs and requeue
those jobs it detects as failed.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.foreman.main import monitor_jobs

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        monitor_jobs()

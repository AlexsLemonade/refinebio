"""
This command will run the Foreman's main function monitor_jobs.
This will cause the Foreman to check for a number of different
failures for both the DownloaderJobs and ProcessorJobs and requeue
those jobs it detects as failed.
"""

from django.core.management.base import BaseCommand
from data_refinery_foreman.foreman.main import monitor_jobs


class Command(BaseCommand):
    def handle(self, *args, **options):
        monitor_jobs()

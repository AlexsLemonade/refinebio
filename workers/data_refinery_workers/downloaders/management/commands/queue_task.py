from django.core.management.base import BaseCommand
from data_refinery_workers.downloaders.array_express \
    import download_array_express


# Just a temporary way to queue a celery task
# without running the surveyor.
class Command(BaseCommand):
    def handle(self, *args, **options):
        download_array_express(0)

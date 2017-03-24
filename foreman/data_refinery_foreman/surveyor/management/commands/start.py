from django.core.management.base import BaseCommand, CommandError
from data_refinery_foreman.surveyor import surveyor


class Command(BaseCommand):
    def handle(self, *args, **options):
        surveyor.test()

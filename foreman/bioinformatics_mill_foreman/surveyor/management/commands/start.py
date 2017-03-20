from django.core.management.base import BaseCommand, CommandError
from bioinformatics_mill_foreman.surveyor import start


class Command(BaseCommand):
    def handle(self, *args, **options):
        start.go()

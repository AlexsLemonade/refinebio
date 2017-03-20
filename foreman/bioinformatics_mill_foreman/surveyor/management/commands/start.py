from django.core.management.base import BaseCommand, CommandError
from bioinformatics_mill_foreman.surveyor import surveyor


class Command(BaseCommand):
    def handle(self, *args, **options):
        surveyor.go()

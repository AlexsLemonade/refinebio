import datetime

from django.core.management.base import BaseCommand
from django_elasticsearch_dsl.registries import registry
from django.utils import timezone


# We'll update for the past 30 minutes every 20 minutes.
UPDATE_WINDOW = datetime.timedelta(minutes=30)


class Command(BaseCommand):
    help = 'Manage elasticsearch index.'

    def handle(self, *args, **options):
        """This command is based off of the 'populate' command of Django ES DSL:

        https://github.com/sabricot/django-elasticsearch-dsl/blob/f6b2e0694e4ed69826c824196ccec5863874c856/django_elasticsearch_dsl/management/commands/search_index.py#L86

        We have updated it so that it will do incremental updates
        rather than looping over the full queryset every time.
        """
        models = set(registry.get_models())

        for doc in registry.get_documents(models):
            start_time = timezone.now() - UPDATE_WINDOW
            qs = doc().get_queryset().filter(last_modified__gt=start_time).order_by('id')
            self.stdout.write("Indexing {} '{}' objects".format(
                qs.count(), qs.model.__name__)
            )
            doc().update(qs)

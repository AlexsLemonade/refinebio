import datetime

from django.core.management.base import BaseCommand
from django.utils import timezone

from django_opensearch_dsl.registries import registry

# We'll update for the past 30 minutes every 20 minutes.
UPDATE_WINDOW = datetime.timedelta(minutes=30)


class Command(BaseCommand):
    help = "Incrementally update the search index for recently-modified docs."

    def handle(self, *args, **options):
        """Incremental counterpart to `opensearch index rebuild`: iterates each
        registered Document and reindexes only the last UPDATE_WINDOW of
        modifications."""
        models = set(registry.get_models())

        for doc in registry.get_documents(models):
            start_time = timezone.now() - UPDATE_WINDOW
            qs = doc().get_queryset().filter(last_modified__gt=start_time).order_by("id")
            self.stdout.write("Indexing {} '{}' objects".format(qs.count(), qs.model.__name__))
            doc().update(qs)

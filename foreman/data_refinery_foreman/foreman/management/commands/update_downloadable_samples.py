from django.core.management.base import BaseCommand
from data_refinery_common.models import Experiment, ComputationalResultAnnotation


class Command(BaseCommand):
    def handle(self, *args, **options):
        organism_ids = list(
            ComputationalResultAnnotation.objects.filter(data__is_qn=True).values_list(
                "data__organism_id", flat=True
            )
        )
        for experiment in Experiment.objects.all():
            experiment.num_downloadable_samples = experiment.samples.filter(
                is_processed=True, organism__id__in=organism_ids
            ).count()
            experiment.save()

        print("Updated the num_downloadable_samples field on all experiment objects.")

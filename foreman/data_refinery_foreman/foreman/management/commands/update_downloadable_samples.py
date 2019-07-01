import django.contrib.postgres.fields
from data_refinery_common.models import Experiment, ComputationalResultAnnotation

class Command(BaseCommand):
    def handle(self, *args, **options):
        for experiment in Experiment.objects.all():
            organism_ids = list(ComputationalResultAnnotation.objects.filter(data__is_qn=True).values_list('data__organism_id', flat=True))
            experiment.num_downloadable_samples = experiment.samples.filter(is_processed=True, organism__id__in=organism_ids).count()
            experiment.save()

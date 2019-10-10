
from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    Sample,
    ComputationalResult,
    ComputationalResultAnnotation,
    Organism
)
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--dry-run',
                            help='Prints resulting actions without actually running them.',
                            action='store_true')
        parser.add_argument("--min",
                            type=int,
                            default=100,
                            help=("Minimum number of processed microarray samples for each organism"))

    def handle(self, *args, **options):
        """ Takes care of removing all qn targets that have less than `min`
        processed microarray samples. """
        remove_invalid_qn_targets(options['min'], options['dry_run'])

def remove_invalid_qn_targets(min_samples, dry_run):
    computational_result_annotations = ComputationalResultAnnotation.objects\
        .filter(data__is_qn=True)\
        .order_by('-created_at')
    qn_target_ids = []
    for annotation in computational_result_annotations:
        organism = Organism.objects.get(id=annotation.data.organism_id)
        samples = Sample.objects.filter(organism=organism, has_raw=True, technology="MICROARRAY", is_processed=True)
        if samples.count() < min_samples:
            # remove the referenced QN Target because it doesn't have enough samples
            logger.debug('Removing QN Target because it did not have enough microarray samples',
                         computational_result=annotation.result,
                         organism=organism)
            qn_target_ids.append(annotation.result.id)

    if not dry_run:
        # delete all invalid qn targets
        ComputationalResult.objects.filter(id__in=qn_target_ids).delete()
    else:
        logger.info("Would have removed computational results with ids %s", str(qn_target_ids))

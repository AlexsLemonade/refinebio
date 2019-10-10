
from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    Sample,
    ComputationalResult,
    ComputationalResultAnnotation,
    Organism
)
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

# Minimum number of microarray samples required to create a QN Target
QN_TARGET_SAMPLE_THRESHOLD = 100

def remove_invalid_qn_targets(dry_run):
    computational_result_annotations = ComputationalResultAnnotation.objects\
        .filter(data__is_qn=True)\
        .order_by('-created_at')
    qn_target_ids = []
    for annotation in computational_result_annotations:
        organism = Organism.objects.get(id=annotation.data.organism_id)
        samples = Sample.objects.filter(organism=organism, has_raw=True, technology="MICROARRAY", is_processed=True)
        if samples.count() < QN_TARGET_SAMPLE_THRESHOLD:
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

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--dry-run',
                            help='Prints resulting actions without actually running them.',
                            action='store_true')
        parser.add_argument('--qn-targets',
                            action='store_true',
                            help='Finds any invalid QN Targets that got generated and removes them from the db.')

    def handle(self, *args, **options):
        """ Re-queues all unprocessed RNA-Seq samples for an organism. """
        if options["qn_targets"]:
            remove_invalid_qn_targets(options['dry_run'])


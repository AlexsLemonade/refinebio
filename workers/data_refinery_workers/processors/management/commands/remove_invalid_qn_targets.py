
from django.core.management.base import BaseCommand
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    Organism
)
from data_refinery_common.logging import get_and_configure_logger
from .create_qn_target import organism_can_have_qn_target
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
    qn_target_ids = []

    for organism in Organism.object.filter(qn_target__isnull=False):
        if not organism_can_have_qn_target(organism, min_samples):
            # Remove all qn targets associated with this object
            qn_target_ids += ComputationalResultAnnotation.objects\
                                .filter(data__is_qn=True, organism_id=organism.id)\
                                .values_list('result__id', flat=True)
            logger.debug('Remove all QN targets for organism', organism=organism)

    if not dry_run:
        # delete all invalid qn targets from S3
        qn_targets = ComputationalResult.objects.filter(id__in=qn_target_ids)
        for qn_target in qn_targets:
            for computed_file in qn_target.computedfile_set.all():
                computed_file.delete_s3_file()

        ComputationalResult.objects.filter(id__in=qn_target_ids).delete()
    else:
        logger.info("Would have removed computational results with ids %s", str(qn_target_ids))
        computed_file_ids = ComputationalResult.objects\
            .filter(id__in=qn_target_ids)\
            .values_list('computedfile_set__id', flat=True)
        logger.info("Would have removed computed files with ids %s", str(computed_file_ids))

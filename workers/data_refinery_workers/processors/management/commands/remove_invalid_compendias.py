
from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (ComputationalResult, ComputedFile,
                                         Organism)

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
        """ Takes care of removing all compendias that are associated with organisms with
        invalid QN targets. """
        remove_invalid_compendias(options['min'], options['dry_run'])


def remove_invalid_compendias(min_samples, dry_run):
    organism_with_compendia = _get_organisms_with_compendias()

    computational_result_ids = []

    for organism in organism_with_compendia:
        if not organism_can_have_qn_target(organism, min_samples):
            # Remove all compendias that are associated with organisms that can't have QN targets
            computational_result_ids += ComputedFile.objects\
                .filter(is_compendia=True, compendia_organism=organism)\
                .values_list('result__id', flat=True)
            logger.debug('Remove all Compendias for organism', organism=organism)

    if not dry_run:
        # delete all invalid compendias from S3
        compendias = ComputationalResult.objects.filter(id__in=computational_result_ids)
        for computational_result in compendias:
            computational_result.remove_computed_files_from_s3()

        # delete all compendias
        compendias.delete()
    else:
        logger.info("Would have removed computational results with ids %s",
                    str(computational_result_ids))


def _get_organisms_with_compendias():
    organism_ids = ComputedFile.objects\
        .filter(is_compendia=True)\
        .values_list('compendia_organism__id', flat=True)\
        .distinct()
    return Organism.objects.filter(id__in=organism_ids)

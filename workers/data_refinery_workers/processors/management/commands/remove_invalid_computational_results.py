from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Organism,
)

from .create_qn_target import organism_can_have_qn_target

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            default=False,
            action="store_true",
            help="Prints resulting actions without actually running them.",
        )
        parser.add_argument(
            "--min",
            type=int,
            default=100,
            help=("Minimum number of processed microarray samples for each organism"),
        )
        parser.add_argument(
            "--qn-target",
            default=False,
            action="store_true",
            help=("Remove invalid QN targets"),
        )
        parser.add_argument(
            "--compendias",
            default=False,
            action="store_true",
            help=("Remove invalid Compendias"),
        )

    def handle(self, *args, **options):
        """ Takes care of removing invalid computational results. """
        computational_result_ids = []

        if options["compendias"]:
            computational_result_ids += remove_invalid_compendias(options["min"])

        if options["qn_target"]:
            computational_result_ids += remove_invalid_qn_targets(options["min"])

        if not computational_result_ids:
            logger.info(
                "Nothing removed. Use options --compendia or --qn-target to select which computational results to check."
            )

        logger.info(
            "Removing computational results with ids %s", str(computational_result_ids)
        )

        if not options["dry_run"]:
            # delete all invalid compendias from S3
            compendias = ComputationalResult.objects.filter(
                id__in=computational_result_ids
            )
            for computational_result in compendias:
                computational_result.remove_computed_files_from_s3()
            # delete all compendias
            compendias.delete()


def remove_invalid_qn_targets(min_samples):
    qn_target_ids = []

    for organism in Organism.object.filter(qn_target__isnull=False):
        if not organism_can_have_qn_target(organism, min_samples):
            # Remove all qn targets associated with this object
            qn_target_ids += ComputationalResultAnnotation.objects.filter(
                data__is_qn=True, data__organism_id=organism.id
            ).values_list("result__id", flat=True)
            logger.debug("Remove all QN targets for organism", organism=organism)

    return qn_target_ids


def remove_invalid_compendias(min_samples):
    organism_with_compendia = _get_organisms_with_compendias()

    computational_result_ids = []

    for organism in organism_with_compendia:
        if not organism_can_have_qn_target(organism, min_samples):
            # Remove all compendias that are associated with organisms that can't have QN targets
            computational_result_ids += ComputedFile.objects.filter(
                is_compendia=True, compendia_organism=organism
            ).values_list("result__id", flat=True)
            logger.debug("Remove all Compendias for organism", organism=organism)

    return computational_result_ids


def _get_organisms_with_compendias():
    organism_ids = (
        ComputedFile.objects.filter(is_compendia=True)
        .values_list("compendia_organism__id", flat=True)
        .distinct()
    )
    return Organism.objects.filter(id__in=organism_ids)

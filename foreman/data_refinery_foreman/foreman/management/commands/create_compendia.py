import sys
from typing import List

from django.core.management.base import BaseCommand
from django.db.models.aggregates import Count
from django.db.models.expressions import Q

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (Dataset, Experiment, Organism,
                                         ProcessorJob,
                                         ProcessorJobDatasetAssociation)
from data_refinery_common.utils import queryset_iterator

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument(
            "--organisms",
            type=str,
            help=("Comma separated list of organism names."))

        parser.add_argument(
            "--quant-sf-only",
            type=lambda x: x == "True",
            help=("Whether to create a quantpendium or normal compendium. Quantpendium will be aggregated by EXPERIMENT"))

        parser.add_argument(
            "--svd-algorithm",
            type=str,
            help=("Specify SVD algorithm applied during imputation ARPACK, RANDOMIZED or NONE to skip."))

    def handle(self, *args, **options):
        """Create a compendium for one or more organisms."""
        target_organisms = self._get_target_organisms(options)
        grouped_organisms = group_organisms_by_biggest_platform(target_organisms)

        quant_sf_only = options["quant_sf_only"] is True

        # default algorithm to arpack until we decide that ranomized is preferred
        svd_algorithm = 'NONE' if quant_sf_only else 'ARPACK'
        if options["svd_algorithm"] in ['ARPACK', 'RANDOMIZED', 'NONE']:
            svd_algorithm = options["svd_algorithm"]

        logger.debug('Generating compendia for organisms', organism_groups=str(grouped_organisms))

        for organism_group in grouped_organisms:
            job = create_job_for_organism(organism_group, quant_sf_only, svd_algorithm)
            logger.info("Sending CREATE_COMPENDIA for Organism", job_id=str(job.pk),
                        organism_group=str(organism_group))
            send_job(ProcessorPipeline.CREATE_COMPENDIA, job)

        sys.exit(0)

    def _get_target_organisms(self, options):
        all_organisms = list(get_compendia_organisms().values_list('id', flat=True))
        all_organisms = Organism.objects.filter(id__in=all_organisms)

        if options["organisms"] is None:
            target_organisms = all_organisms.exclude(name__in=["HOMO_SAPIENS", "MUS_MUSCULUS"])
        else:
            organisms = options["organisms"].upper().replace(" ", "_").split(",")
            target_organisms = all_organisms.filter(name__in=organisms)

        return target_organisms


def get_compendia_organisms():
    """ We start with the organisms that have QN targets associated with them.
    Here we are filtering out the ones with invalid QN targets. This can be removed once
    https://github.com/AlexsLemonade/refinebio/issues/1757 is addressed """
    return Organism.objects\
        .filter(qn_target__isnull=False)\
        .annotate(num_samples=Count('sample', filter=Q(sample__is_processed=True, sample__technology='MICROARRAY'), distinct=True))\
        .filter(num_samples__gt=100)


def group_organisms_by_biggest_platform(target_organisms):
    """ Create groups of organisms that share the same platform
    ref https://github.com/AlexsLemonade/refinebio/issues/1736 """
    result = []

    # process the organisms with the most platforms first
    target_organisms_sorted = target_organisms\
        .annotate(num_platforms=Count('sample__platform_accession_code', distinct=True))\
        .order_by('-num_platforms')

    for organism in target_organisms_sorted:
        added_with_another_organism = any(
            any(item == organism for item in group) for group in result)
        if added_with_another_organism:
            continue

        organism_group = [organism]

        platform_list = list(get_organism_microarray_platforms(organism))

        genus_prefix = organism.get_genus() + '_'
        organisms_sharing_genus = Organism.objects.filter(name__startswith=genus_prefix)
        for genus_organism in organisms_sharing_genus:
            # Group together the orgenisms that share at least one microarray platform
            platforms = get_organism_microarray_platforms(genus_organism)\
                .filter(platform_accession_code__in=platform_list)
            if platforms.exists():
                organism_group.append(genus_organism)

        result.append(organism_group)

    return result


def get_organism_microarray_platforms(organism):
    """ Returns the accession codes of the Affymetrix microarray platforms associated with an
    organism. Ordered by the number of samples for each platform in descending order """
    return organism.sample_set.filter(has_raw=True, technology="MICROARRAY", is_processed=True)\
                   .values('platform_accession_code')\
                   .annotate(count=Count('id'))\
                   .order_by('-count')\
                   .values_list('platform_accession_code', flat=True)


def create_job_for_organism(organisms: List[Organism], quant_sf_only=False, svd_algorithm='ARPACK'):
    """Returns a compendia job for the provided organism.

    Fetch all of the experiments and compile large but normally formated Dataset.
    """
    job = ProcessorJob()
    job.pipeline_applied = ProcessorPipeline.CREATE_COMPENDIA.value
    job.save()

    dataset = Dataset()
    dataset.data = get_dataset(organisms)
    dataset.scale_by = 'NONE'
    # The quantpendias should be aggregated by species
    dataset.aggregate_by = 'EXPERIMENT' if quant_sf_only else 'SPECIES'
    dataset.quantile_normalize = False
    dataset.quant_sf_only = quant_sf_only
    dataset.svd_algorithm = svd_algorithm
    dataset.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = job
    pjda.dataset = dataset
    pjda.save()

    return job


def get_dataset(organisms: List[Organism]):
    """ Builds a dataset with the samples associated with the given organisms """
    dataset = {}

    filter_query = Q()
    for organism in organisms:
        filter_query = filter_query | Q(organisms=organism)

    experiments = Experiment.objects.filter(filter_query).prefetch_related('samples')

    for experiment in queryset_iterator(experiments):
        filter_query = Q()
        for organism in organisms:
            filter_query = filter_query | Q(organism=organism)

        dataset[experiment.accession_code] = list(experiment.samples.filter(filter_query)
                                                  .values_list('accession_code', flat=True))

    return dataset

import sys

from django.core.management.base import BaseCommand

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (Dataset, Experiment, Organism,
                                         ProcessorJob,
                                         ProcessorJobDatasetAssociation)
from data_refinery_common.utils import queryset_iterator

logger = get_and_configure_logger(__name__)

def create_job_for_organism(organism=Organism, quant_sf_only=False, svd_algorithm='ARPACK'):
    """Returns a compendia job for the provided organism.

    Fetch all of the experiments and compile large but normally formated Dataset.
    """
    data = {}
    experiments = Experiment.objects.filter(organisms=organism).prefetch_related('samples')

    for experiment in queryset_iterator(experiments):
        data[experiment.accession_code] = list(experiment.samples.filter(organism=organism)\
            .values_list('accession_code', flat=True))

    job = ProcessorJob()
    if quant_sf_only:
        job.pipeline_applied = ProcessorPipeline.CREATE_QUANTPENDIA.value
    else:
        job.pipeline_applied = ProcessorPipeline.CREATE_COMPENDIA.value
    job.save()

    dset = Dataset()
    dset.data = data
    dset.scale_by = 'NONE'
    dset.aggregate_by = 'EXPERIMENT' if quant_sf_only else 'SPECIES'
    dset.quantile_normalize = False
    dset.quant_sf_only = quant_sf_only
    dset.svd_algorithm = svd_algorithm
    dset.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = job
    pjda.dataset = dset
    pjda.save()

    return job


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
        """Create a compendium for one or more organisms.

        If --organism is supplied will immediately create a compedium
        for it. If not a new job will be dispatched for each organism
        with enough microarray samples except for human and mouse.
        """
        if options["organisms"] is None:
            all_organisms = Organism.objects.exclude(name__in=["HOMO_SAPIENS", "MUS_MUSCULUS"])
        else:
            organisms = options["organisms"].upper().replace(" ", "_").split(",")
            all_organisms = Organism.objects.filter(name__in=organisms)

        # I think we could just use options["quant_sf_only"] but I
        # wanna make sure that values that are not True do not trigger
        # a truthy evaluation.
        quant_sf_only = options["quant_sf_only"] is True

        # default algorithm to arpack until we decide that ranomized is preferred
        svd_algorithm = 'NONE' if quant_sf_only else 'ARPACK'
        if options["svd_algorithm"] in ['ARPACK', 'RANDOMIZED', 'NONE']:
            svd_algorithm = options["svd_algorithm"]

        logger.debug('Generating compendia for organisms', organisms=all_organisms)

        for organism in all_organisms:
            job = create_job_for_organism(organism, quant_sf_only, svd_algorithm)
            logger.info("Sending CREATE_COMPENDIA for Organism", job_id=str(job.pk), organism=str(organism))
            send_job(ProcessorPipeline.CREATE_COMPENDIA, job)

        sys.exit(0)

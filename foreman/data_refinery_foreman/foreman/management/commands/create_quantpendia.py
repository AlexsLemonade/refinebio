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

def create_job_for_organism(organism: Organism):
    """Returns a quantpendia job for the provided organism."""
    data = {}
    experiments = Experiment.objects.filter(
            organisms=organism,
            samples__results__computedfile__filename='quant.sf'
        )\
        .distinct()\
        .prefetch_related('samples')

    for experiment in queryset_iterator(experiments):
        # only include the samples from the target organism that have quant.sf files
        samples_with_quantsf = experiment.samples\
            .filter(
                organism=organism,
                results__computedfile__filename='quant.sf'
            )\
            .values_list('accession_code', flat=True)\
            .distinct()
        data[experiment.accession_code] = list(samples_with_quantsf)

    job = ProcessorJob()
    job.pipeline_applied = ProcessorPipeline.CREATE_QUANTPENDIA.value
    job.save()

    dset = Dataset()
    dset.data = data
    dset.scale_by = 'NONE'
    dset.aggregate_by = 'EXPERIMENT'
    dset.quantile_normalize = False
    dset.quant_sf_only = True
    dset.svd_algorithm = 'NONE'
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

    def handle(self, *args, **options):
        """Create a quantpendia for one or more organisms."""
        all_organisms = Organism.objects.all().filter(qn_target__isnull=False)
        if options["organisms"] is None:
            all_organisms = all_organisms.exclude(name__in=["HOMO_SAPIENS", "MUS_MUSCULUS"])
        else:
            organisms = options["organisms"].upper().replace(" ", "_").split(",")
            all_organisms = all_organisms.filter(name__in=organisms)

        logger.debug('Generating quantpendia for organisms', organisms=all_organisms)

        for organism in all_organisms:
            job = create_job_for_organism(organism)
            logger.info("Sending compendia job for Organism", job_id=str(job.pk), organism=str(organism))
            send_job(ProcessorPipeline.CREATE_QUANTPENDIA, job)

        sys.exit(0)

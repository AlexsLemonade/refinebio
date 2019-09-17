import sys

from django.core.management.base import BaseCommand
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job

from data_refinery_common.models import (
    Experiment,
    Sample,
    Organism,
    ProcessorJob,
    Dataset,
    ProcessorJobDatasetAssociation,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation
)
from data_refinery_workers.processors import create_compendia

logger = get_and_configure_logger(__name__)


def create_job_for_organism(organism=Organism):
    """Returns a compendia job for the provided organism.

    Fetch all of the experiments and compile large but normally formated Dataset.
    """
    data = {}
    experiments = Experiment.objects.filter(organisms=organism).prefetch_related('samples')
    for experiment in experiments:
        data[experiment.accession_code] = list(experiment.samples.filter(organism=organism).values_list('accession_code', flat=True))

    job = ProcessorJob()
    job.pipeline_applied = "COMPENDIA"
    job.save()

    dset = Dataset()
    dset.data = data
    dset.scale_by = 'NONE'
    dset.aggregate_by = 'SPECIES'
    dset.quantile_normalize = False
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
        """Create a compendium for one or more organisms.

        If --organism is supplied will immediately create a compedium
        for it. If not a new job will be dispatched for each organism
        with enough microarray samples.
        """
        if options["organisms"] is None:
            all_organisms = Organism.objects.all()
        else:
            organisms = options["organisms"].upper().replace(" ", "_").split(",")
            all_organisms = Organism.objects.filter(name__in=organisms)

        logger.error(all_organisms)

        if all_organisms.count() > 1:
            for organism in all_organisms:
                logger.error(organism)
                job = create_job_for_organism(organism)
                logger.info("Sending CREATE_COMPENDIA for Organism", job_id=str(job.pk), organism=str(organism))
                send_job(ProcessorPipeline.CREATE_COMPENDIA, job)
        else:
            job = create_job_for_organism(all_organisms[0])
            create_compendia.create_compendia(job.id)

        sys.exit(0)

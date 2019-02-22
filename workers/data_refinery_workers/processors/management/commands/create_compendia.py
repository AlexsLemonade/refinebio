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
    OrganismIndex,
    ExperimentSampleAssociation
)
from data_refinery_workers.processors import create_compendia

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument(
            "--organism",
            type=str,
            help=("Name of organism"))

    def handle(self, *args, **options):
        """ For every (or a supplied) organism, fetch all of the experiments and compile large but normally formated Dataset.

        Send all of them to the Smasher. Smash them. Retrieve manually as desired.
        """

        dataset_ids = []

        if options["organism"] is None:
            all_organisms = Organism.objects.all()
        else:
            all_organisms = [Organism.get_object_for_name(options["organism"].upper())]

        for organism in all_organisms:
            data = {}
            experiments = Experiment.objects.filter(id__in=(ExperimentOrganismAssociation.objects.filter(organism=organism)).values('experiment'))
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

            final_context = create_compendia.create_compendia(job.id)

        sys.exit(0)

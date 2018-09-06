import requests
import sys

from django.core.management.base import BaseCommand
from django.db.models import Count

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SampleComputedFileAssociation,
)
from data_refinery_workers.processors import qn_reference, utils


logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--organism",
            type=str,
            help=("Name of organism"))

    def handle(self, *args, **options):
        """
        """

        if options["organism"] is None:
            logger.error("You must specify an organism")
            sys.exit(1)

        organism = Organism.get_object_for_name(options["organism"].upper())
        samples = Sample.processed_objects.filter(organism=organism)
        if samples.count() == 0:
            logger.error("No processed samples for that organism.")
            sys.exit(1)

        platform_counts = samples.values('platform_accession_code').annotate(dcount=Count('platform_accession_code')).order_by('-dcount')

        biggest_platform = platform_counts[0]['platform_accession_code']
        sample_codes_results = Sample.objects.filter(platform_accession_code=biggest_platform).values('accession_code')
        sample_codes = [res['accession_code'] for res in sample_codes_results]

        dataset = Dataset()
        dataset.data = {organism.name + '_(' + biggest_platform + ')': sample_codes}
        dataset.aggregate_by = "ALL"
        dataset.scale_by = "NONE"
        dataset.quantile_normalize = False
        dataset.save()

        job = ProcessorJob()
        job.pipeline_applied = "QN_REFERENCE"
        job.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        final_context = qn_reference.create_qn_reference(job.pk)

        if final_context['success']:
            print(":D")
            sys.exit(0)
        else:
            print(":(")
            sys.exit(-1)

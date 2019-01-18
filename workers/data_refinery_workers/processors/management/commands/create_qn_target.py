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
        parser.add_argument(
            "--all",
            default=False,
            action='store_true',
            help=("Perform for all organisms that meet <min> processed samples"))
        parser.add_argument(
            "--min",
            type=int,
            default=100,
            help=("Minimum number of processed samples"))
        parser.add_argument(
            "--platform",
            type=str,
            default=None,
            help=("Name of platform")
            )

    def handle(self, *args, **options):
        """
        """

        if options["organism"] is None and not options["all"]:
            logger.error("You must specify an organism or --all")
            sys.exit(1)

        if options["organism"]:
            organisms = [Organism.get_object_for_name(options["organism"].upper())]
        else:
            organisms = Organism.objects.all()

        for organism in organisms:
            samples = Sample.processed_objects.filter(organism=organism, has_raw=True, technology="MICROARRAY", is_processed=True)
            if samples.count() == 0:
                logger.error("No processed samples for organism.",
                    organism=organism,
                    count=samples.count()
                    )
                continue
            if samples.count() < options['min']:
                logger.error("Proccessed samples don't meet minimum threshhold",
                    organism=organism,
                    count=samples.count(),
                    min=options["min"]
                )
                continue

            if options["platform"] is None:
                platform_counts = samples.values('platform_accession_code').annotate(dcount=Count('platform_accession_code')).order_by('-dcount')
                biggest_platform = platform_counts[0]['platform_accession_code']
            else:
                biggest_platform = options["platform"]

            sample_codes_results = Sample.processed_objects.filter(
                platform_accession_code=biggest_platform, 
                has_raw=True, 
                technology="MICROARRAY", 
                is_processed=True).values('accession_code')
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
                self.stdout.write("Target file: " + final_context['target_file'])
                self.stdout.write("Target S3: " + str(final_context['computed_files'][0].get_s3_url()))
            else:
                print(":(")

import sys
import requests

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

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--dataset-ids",
            type=str,
            help=("Comma-separated Dataset IDs."))

    def handle(self, *args, **options):
        """ For every organism, fetch all of the experiments and compile large but normally formated Dataset.

        Send all of them to the Smasher. Smash them. Retrieve manually as desired.
        """

        if options["dataset_ids"] is None:
            logger.error("You must specify dataset IDs.")
            sys.exit(1) 

        dataset_ids = options["dataset_ids"].split(',')

        datasets = Dataset.objects.filter(id__in=dataset_ids)
        for dataset in datasets:
            try:
                organism = dataset.email_address.split('@')[0].split('compendia_')[1]
                organism_name = organism
            except:
                organism_name = str(dataset.pk)

            if dataset.is_available:
                print(organism_name + ": " + dataset.s3_url())
            
                response = requests.get(dataset.s3_url(), stream=True)
                response.raise_for_status()
                with open(organism_name + '.zip', 'wb') as handle:
                    for block in response.iter_content(1024):
                        handle.write(block)
            else:
                print(organism_name + " is not available!")

        sys.exit(0)

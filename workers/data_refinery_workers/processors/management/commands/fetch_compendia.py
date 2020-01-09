import sys

from django.core.management.base import BaseCommand

import requests

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    Dataset,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
)

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--dataset-ids", type=str, help=("Comma-separated Dataset IDs."))

    def handle(self, *args, **options):
        """ Given a dataset ID, fetch the resulting smashed object.
        """

        if options["dataset_ids"] is None:
            logger.error("You must specify dataset IDs.")
            sys.exit(1)

        dataset_ids = options["dataset_ids"].split(",")

        datasets = Dataset.objects.filter(id__in=dataset_ids)
        for dataset in datasets:
            try:
                organism_name = dataset.email_address.split("@")[0].split("compendia_")[1]
            except:
                organism_name = str(dataset.pk)

            if dataset.is_available:
                print(organism_name + ": " + dataset.s3_url())

                response = requests.get(dataset.s3_url(), stream=True)
                response.raise_for_status()
                with open(organism_name + ".zip", "wb") as handle:
                    for block in response.iter_content(1024):
                        handle.write(block)
            else:
                print(organism_name + " is not available!")

        sys.exit(0)

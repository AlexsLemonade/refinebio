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

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):

    def handle(self, *args, **options):
        """ For every organism, fetch all of the experiments and compile large but normally formated Dataset.

        Send all of them to the Smasher. Smash them. Retrieve manually as desired.
        """

        dataset_ids = []

        all_organisms = Organism.objects.all()
        for organism in all_organisms:
            data = {}
            experiments = Experiment.objects.filter(id__in=(ExperimentOrganismAssociation.objects.filter(organism=organism)).values('experiment'))
            for experiment in experiments:
                data[experiment.accession_code] = list(experiment.samples.values_list('accession_code', flat=True))

            dataset = Dataset()
            dataset.data = data
            dataset.aggregate_by = "EXPERIMENT"
            dataset.scale_by = "MINMAX"
            dataset.email_address = "ccdl_compendia_" + organism.name + "@mailinator.com"
            dataset.save()
            dataset_ids.append(str(dataset.pk))

            processor_job = ProcessorJob()
            processor_job.pipeline_applied = "SMASHER"
            processor_job.ram_amount = 4096
            processor_job.save()

            pjda = ProcessorJobDatasetAssociation()
            pjda.processor_job = processor_job
            pjda.dataset = dataset
            pjda.save()

            try:
                send_job(ProcessorPipeline.SMASHER, processor_job)
            except Exception as e:
                print("Couldn't send job for organism: " + organism.name)
                print(str(e))
                continue

            print("Created compendia smasher job for " + organism.name + " with " + str(len(data.keys())) + " experiments!")

        print("Smashers dispatched, to retrieve:")
        print("python manage.py fetch_compendia --dataset-ids " + ','.join(dataset_ids))

        sys.exit(0)

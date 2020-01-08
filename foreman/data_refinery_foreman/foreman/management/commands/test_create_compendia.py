from django.core.management import call_command
from django.test import TransactionTestCase

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)

from .test_create_quantpendia import get_organism_with_qn_target, make_test_data


class CompendiaCommandTestCase(TransactionTestCase):
    def test_compendia_command(self):
        organism = get_organism_with_qn_target()
        make_test_data(organism)

        try:
            call_command("create_compendia", organisms=organism.name)
        except SystemExit as e:  # this is okay!
            pass

        processor_job = (
            ProcessorJob.objects.filter(pipeline_applied="CREATE_COMPENDIA")
            .order_by("-created_at")
            .first()
        )

        # check that the processor job was created correctly
        self.assertIsNotNone(processor_job)
        self.assertEquals(processor_job.datasets.first().data, {"GSE51088": ["GSM1237818"]})

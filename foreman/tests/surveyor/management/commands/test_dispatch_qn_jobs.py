from unittest.mock import patch

from django.test import TransactionTestCase, tag

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
)
from data_refinery_common.models.organism import Organism
from data_refinery_foreman.surveyor.management.commands.dispatch_qn_jobs import Command


class QNRefTestCase(TransactionTestCase):
    @tag("foreman")
    @patch("data_refinery_foreman.surveyor.management.commands.dispatch_qn_jobs.send_job")
    def test_qn_reference(self, mock_send_job):
        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        organism.save()

        experiment = Experiment()
        experiment.accession_code = "12345"
        experiment.save()

        for code in [str(i) for i in range(1, 401)]:
            sample = Sample()
            sample.accession_code = code
            sample.title = code
            sample.platform_name = f"Affymetrix {organism.name}"
            sample.platform_accession_code = f"A-MEXP-{organism.name}"
            sample.manufacturer = "AFFYMETRIX"
            sample.organism = organism
            sample.technology = "MICROARRAY"
            sample.is_processed = True
            sample.has_raw = True
            sample.save()

            cr = ComputationalResult()
            cr.save()

            computed_file = ComputedFile()
            computed_file.filename = code + ".tsv"
            computed_file.absolute_file_path = "/home/user/data_store/QN/" + code + ".tsv"
            computed_file.size_in_bytes = int(code)
            computed_file.result = cr
            computed_file.is_smashable = True
            computed_file.save()

            scfa = SampleComputedFileAssociation()
            scfa.sample = sample
            scfa.computed_file = computed_file
            scfa.save()

            exsa = ExperimentSampleAssociation()
            exsa.experiment = experiment
            exsa.sample = sample
            exsa.save()

            # We need more than one organism for the tests, but can't
            # repeat accesion codes, so halfway through just change the organism.
            if int(code) == 200:
                organism = Organism(name="MUS_MUSCULUS", taxonomy_id=111)
                organism.save()

        # Setup is done, actually run the command.
        command = Command()
        command.handle(organisms="HOMO_SAPIENS,MUS_MUSCULUS")

        self.assertEqual(len(mock_send_job.mock_calls), 2)
        self.assertEqual(ProcessorJob.objects.count(), 2)

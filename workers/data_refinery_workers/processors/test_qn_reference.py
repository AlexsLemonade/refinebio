import os

from django.core.management import call_command
from django.test import TransactionTestCase, tag
from django.utils.six import StringIO

from data_refinery_common.models import (
    ProcessorJob,
    Sample,
    SampleComputedFileAssociation,
    Dataset,
    ComputedFile,
    ComputationalResult,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJobDatasetAssociation
)
from data_refinery_common.models.organism import Organism
from data_refinery_workers.processors import qn_reference, smasher, utils


class QNRefTestCase(TransactionTestCase):
    @tag('qn')
    def test_qn_reference(self):
        job = ProcessorJob()
        job.pipeline_applied = "QN_REFERENCE"
        job.save()

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        homo_sapiens.save()

        experiment = Experiment()
        experiment.accession_code = "12345"
        experiment.organism_names = [homo_sapiens.name]
        experiment.save()
        # We don't have a 0.tsv
        codes = [str(i) for i in range(1, 201)]

        for code in codes:
            sample = Sample()
            sample.accession_code = code
            sample.title = code
            sample.platform_accession_code = 'A-MEXP-1171'
            sample.manufacturer = "SLIPPERY DICK'S DISCOUNT MICROARRAYS"
            sample.organism = homo_sapiens
            sample.technology = "MICROARRAY"
            sample.is_processed = True
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

        dataset = Dataset()
        dataset.data = {"12345": ["1", "2", "3", "4", "5", "6"]}
        dataset.aggregate_by = "ALL"
        dataset.scale_by = "NONE"
        dataset.quantile_normalize = False  # We don't QN because we're creating the target now
        dataset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        final_context = qn_reference.create_qn_reference(job.pk)
        self.assertTrue(final_context['success'])
        self.assertTrue(os.path.exists(final_context['target_file']))
        self.assertEqual(os.path.getsize(final_context['target_file']), 562)

        homo_sapiens.refresh_from_db()
        target = homo_sapiens.qn_target.computedfile_set.latest()
        self.assertEqual(target.sha1, '636d72d5cbf4b9785b0bd271a1430b615feaa7ea')

        # Create and run a smasher job that will use the QN target we just made.
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        ds = Dataset()
        ds.data = {"12345": ["1", "2", "3", "4", "5"]}
        ds.aggregate_by = 'SPECIES'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = True
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(final_context['success'])

        self.assertEqual(final_context['merged_qn']['1'][0], -0.4379488528812934)
        self.assertEqual(final_context['original_merged']['1'][0], -0.576210936113982)

    @tag('qn')
    def test_qn_management_command(self):
        """Test that the management command fires off and then does not create
        a job for an organism that does not have enough samples on the same
        platform."""

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        homo_sapiens.save()

        experiment = Experiment()
        experiment.accession_code = "12345"
        experiment.organism_names = [homo_sapiens.name]
        experiment.save()
        codes = ["1", "2", "3", "4", "5", "6"]
        # We don't have a 0.tsv

        for code in codes:
            sample = Sample()
            sample.accession_code = code
            sample.title = code
            sample.platform_accession_code = 'A-MEXP-1171'
            sample.manufacturer = "SLIPPERY DICK'S DISCOUNT MICROARRAYS"
            sample.organism = homo_sapiens
            sample.technology = "MICROARRAY"
            sample.is_processed = True
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

        out = StringIO()
        try:
            call_command('create_qn_target', organism='homo_sapiens', min=1, stdout=out)
        except SystemExit as e:  # this is okay!
            pass

        stdout = out.getvalue()
        self.assertFalse('Target file' in stdout)

        # There's not enough samples available in this scenario so we
        # shouldn't have even made a processor job.
        self.assertEqual(ProcessorJob.objects.count(), 0)

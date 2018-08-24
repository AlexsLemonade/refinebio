import os
import shutil
from contextlib import closing
from django.test import TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    Sample,
    Organism,
    SampleComputedFileAssociation,
    ProcessorJobOriginalFileAssociation,
    Dataset,
    ComputedFile,
    ComputationalResult,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJobDatasetAssociation
)
from data_refinery_workers.processors import qn_reference, smasher, utils


class QNRefTestCase(TestCase):

    @tag('qn')
    def test_sanity(self):
        print("Hey!")

    @tag('qn')
    def test_qn_reference(self):
        job = ProcessorJob()
        job.pipeline_applied = "QN_REFERENCE"
        job.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        experiment = Experiment()
        experiment.accession_code = "12345"
        experiment.save()

        for code in ['1', '2', '3', '4', '5']:
            sample = Sample()
            sample.accession_code = code
            sample.title = code
            sample.platform_accession_code = 'A-MEXP-1171'
            sample.manufacturer = "ILLUMINA"
            sample.organism = homo_sapiens
            sample.is_processed = True
            sample.save()

            cr = ComputationalResult()
            cr.save()

            file = ComputedFile()
            file.filename = code + ".tsv"
            file.absolute_file_path = "/home/user/data_store/QN/" + code + ".tsv"
            file.size_in_bytes = int(code)
            file.result = cr
            file.is_smashable = True
            file.save()

            scfa = SampleComputedFileAssociation()
            scfa.sample = sample
            scfa.computed_file = file
            scfa.save()

            exsa = ExperimentSampleAssociation()
            exsa.experiment = experiment
            exsa.sample = sample
            exsa.save()

        
        dataset = Dataset()
        dataset.data = {"12345": ["1", "2", "3", "4", "5"]}
        dataset.aggregate_by = "ALL"
        dataset.scale_by = "NONE"
        dataset.quantile_normalize = False # We don't QN because we're creating the target now
        dataset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        final_context = qn_reference.create_qn_reference(job.pk)
        self.assertTrue(final_context['success'])

        self.assertTrue(os.path.exists(final_context['target_file']))
        self.assertEqual(os.path.getsize(final_context['target_file']), 519)

        target = utils.get_most_recent_qn_target_for_organism(homo_sapiens)
        self.assertEqual(target.sha1, 'a38ae13de860e47e0251dd02d1a8e88f576d83ad')

        ###
        # Smasher with QN
        ###

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

        self.assertEqual(final_context['merged_qn']['1'][0], -0.437948852881293)
        self.assertEqual(final_context['original_merged']['1'][0], -0.576210936113982)

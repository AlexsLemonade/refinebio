import os
import shutil
from contextlib import closing
from django.test import TransactionTestCase, TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.job_lookup import ProcessorPipeline
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
    ComputationalResultAnnotation,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJobDatasetAssociation,
    SampleAnnotation,
    SampleResultAssociation
)
from data_refinery_workers.processors.create_quantpendia import create_quantpendia


class QuantpendiaTestCase(TransactionTestCase):
    @tag("compendia")
    def test_create_quantpendia(self):
        job = ProcessorJob()
        job.pipeline_applied = ProcessorPipeline.CREATE_QUANTPENDIA.value
        job.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51088"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        sample = Sample()
        sample.accession_code = 'GSM1237818'
        sample.title = 'GSM1237818'
        sample.organism = homo_sapiens
        sample.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = experiment
        esa.sample = sample
        esa.save()

        computed_file = ComputedFile()
        computed_file.s3_key = "smasher-test-quant.sf"
        computed_file.s3_bucket = "data-refinery-test-assets"
        computed_file.filename = "quant.sf"
        computed_file.absolute_file_path = "/home/user/data_store/QUANT/smasher-test-quant.sf"
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.size_in_bytes = 123123
        computed_file.sha1 = "08c7ea90b66b52f7cd9d9a569717a1f5f3874967" # this matches with the downloaded file
        computed_file.save()

        computed_file = ComputedFile()
        computed_file.filename = "logquant.tsv"
        computed_file.is_smashable = True
        computed_file.size_in_bytes = 123123
        computed_file.result = result
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {'GSE51088': ['GSM1237818']}
        ds.aggregate_by = 'EXPERIMENT'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.quant_sf_only = True # Make the dataset include quant.sf files only
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = create_quantpendia(job.id)

        self.assertTrue(os.path.exists(final_context['output_dir'] + '/GSE51088/GSM1237818_quant.sf'))
        self.assertTrue(os.path.exists(final_context['output_dir'] + '/README.md'))
        self.assertTrue(os.path.exists(final_context['output_dir'] + '/LICENSE.TXT'))
        self.assertTrue(os.path.exists(final_context['output_dir'] + '/aggregated_metadata.json'))

        self.assertTrue(final_context['metadata']['quant_sf_only'])
        self.assertEqual(final_context['metadata']['num_samples'], 1)
        self.assertEqual(final_context['metadata']['num_experiments'], 1)


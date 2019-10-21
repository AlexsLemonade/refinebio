from django.test import TransactionTestCase
from django.utils import timezone
from django.core.management import call_command

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.models import (ComputationalResult,
                                         ComputedFile,
                                         Dataset,
                                         Experiment,
                                         ExperimentSampleAssociation,
                                         Organism,
                                         ProcessorJob,
                                         ProcessorJobDatasetAssociation,
                                         Sample,
                                         SampleComputedFileAssociation,
                                         SampleResultAssociation)
from data_refinery_foreman.surveyor.test_end_to_end import wait_for_job

from .create_compendia import create_job_for_organism


class CompendiaCommandTestCase(TransactionTestCase):
    def test_quantpendia_command(self):
        self.make_test_data()

        args = []
        options = {'organisms': 'HOMO_SAPIENS', 'quant_sf_only': True}
        call_command('create_compendia', *args, **options)

        start_time = timezone.now()
        processor_job = ProcessorJob.objects.filter(pipeline_applied='CREATE_QUANTPENDIA').first()
        wait_for_job(processor_job, ProcessorJob, start_time)
        self.assertTrue(processor_job.success)

        # assert we have a quantpendia for human
        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        quantpendias = ComputedFile.objects.filter(is_compendia=True, quant_sf_only=True, compendia_organism=homo_sapiens).first()
        self.assertTrue(quantpendias)

    def make_test_data(self):
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

import os
import zipfile

from django.test import TestCase, tag
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    ComputationalResult,
    ComputedFile,
    Experiment,
    Organism,
    Sample,
    SampleResultAssociation,
    ExperimentSampleAssociation,
    Dataset,
    ProcessorJobDatasetAssociation,
    SampleComputedFileAssociation
)
from data_refinery_workers.processors import smasher
from data_refinery_workers.processors import utils

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SMASHER"
    pj.save()

    experiment = Experiment()
    experiment.accession_code = "GSE51081"
    experiment.save()

    result = ComputationalResult()
    result.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

    sample = Sample()
    sample.accession_code = 'GSM1237810'
    sample.title = 'GSM1237810'
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
    computed_file.filename = "GSM1237810_T09-1084.PCL"
    computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.is_smashable = True
    computed_file.save()

    assoc = SampleComputedFileAssociation()
    assoc.sample = sample
    assoc.computed_file = computed_file
    assoc.save()

    sample = Sample()
    sample.accession_code = 'GSM1237812'
    sample.title = 'GSM1237812'
    sample.organism = homo_sapiens
    sample.save()

    esa = ExperimentSampleAssociation()
    esa.experiment = experiment
    esa.sample = sample
    esa.save()

    assoc = SampleComputedFileAssociation()
    assoc.sample = sample
    assoc.computed_file = computed_file
    assoc.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

    computed_file = ComputedFile()
    computed_file.filename = "GSM1237812_S97-PURE.PCL"
    computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.is_smashable = True
    computed_file.save()

    assoc = SampleComputedFileAssociation()
    assoc.sample = sample
    assoc.computed_file = computed_file
    assoc.save()

    ds = Dataset()
    ds.data = {'GSE51081': ['GSM1237810', 'GSM1237812']}
    ds.aggregate_by = 'EXPERIMENT' # [ALL or SPECIES or EXPERIMENT]
    ds.scale_by = 'STANDARD' # [NONE or MINMAX or STANDARD or ROBUST]
    ds.email_address = "null@derp.com"
    #ds.email_address = "miserlou+heyo@gmail.com"
    ds.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = pj
    pjda.dataset = ds
    pjda.save()

    return pj

class SmasherTestCase(TestCase):

    @tag("smasher")
    def test_smasher(self):
        """ Main tester. """
        job = prepare_job()

        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)
        dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
        job_context_check = {}
        job_context_check["dataset"] = dataset
        job_context_check["samples"] = dataset.get_samples()
        job_context_check["experiments"] = dataset.get_experiments()
        self.assertEqual(len(job_context_check['samples']), 2)
        self.assertEqual(len(job_context_check['experiments']), 1)

        # Smoke test while we're here..
        dataset.get_samples_by_experiment()
        dataset.get_samples_by_species()
        dataset.get_aggregated_samples()
        dsid = dataset.id

        # XXX: agg_type 'SPECIES' hangs on Linux, not OSX.
        # Don't know why yet.
        # for ag_type in ['ALL', 'EXPERIMENT', 'SPECIES']:
        for ag_type in ['ALL', 'EXPERIMENT']:
            dataset = Dataset.objects.get(id=dsid)
            dataset.aggregate_by = ag_type
            dataset.save()

            final_context = smasher.smash(job.pk, upload=False)

            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context['output_file']), 0)
            self.assertEqual(final_context['dataset'].is_processed, True)

            dataset = Dataset.objects.get(id=dsid)
            dataset.is_processed = False
            dataset.save()

            # Cleanup
            os.remove(final_context['output_file'])

        for scale_type in ['NONE', 'MINMAX', 'STANDARD', 'ROBUST']:
            dataset = Dataset.objects.get(id=dsid)
            dataset.aggregate_by = 'EXPERIMENT'
            dataset.scale_by = scale_type
            dataset.save()

            print ("Smashing " + scale_type)
            final_context = smasher.smash(job.pk, upload=False)

            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context['output_file']), 0)
            self.assertEqual(final_context['dataset'].is_processed, True)

            dataset = Dataset.objects.get(id=dsid)
            dataset.is_processed = False
            dataset.save()

            # Cleanup
            os.remove(final_context['output_file'])

        # Make sure we can use our outputs mathematically
        for scale_type in ['MINMAX', 'STANDARD']:
            dataset = Dataset.objects.get(id=dsid)
            dataset.aggregate_by = 'EXPERIMENT'
            dataset.scale_by = scale_type
            dataset.save()

            final_context = smasher.smash(job.pk, upload=False)
            final_frame = final_context['final_frame']

            # Sanity test that these frames can be computed upon
            final_frame.mean(axis=1)
            final_frame.min(axis=1)
            final_frame.max(axis=1)
            final_frame.std(axis=1)
            final_frame.median(axis=1)

            zf = zipfile.ZipFile(final_context['output_file'])
            namelist = zf.namelist()

            self.assertTrue('metadata.tsv' in namelist)
            self.assertTrue('metadata.json' in namelist)
            self.assertTrue('README.md' in namelist)
            self.assertTrue('LICENSE.TXT' in namelist)
            self.assertTrue('GSE51081.tsv' in namelist)

            os.remove(final_context['output_file'])


    @tag("smasher")
    def test_get_results(self):
        """ Test our ability to collect the appropriate samples. """

        sample = Sample()
        sample.accession_code = 'GSM45588'
        sample.save()

        result = ComputationalResult()
        result.save()

        computed_file1 = ComputedFile()
        computed_file1.filename = "oh_boy.txt"
        computed_file1.result = result
        computed_file1.size_in_bytes = 123
        computed_file1.is_smashable = True
        computed_file1.save()

        computed_file2 = ComputedFile()
        computed_file2.filename = "gee_whiz.bmp"
        computed_file2.result = result
        computed_file2.size_in_bytes = 123
        computed_file2.is_smashable = False
        computed_file2.save()

        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file1
        assoc.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file2
        assoc.save()

        computed_files = sample.get_result_files()
        self.assertEqual(computed_files.count(), 2)

        smashable_file = sample.get_most_recent_smashable_result_file()
        self.assertTrue(smashable_file.is_smashable)

    @tag("smasher")
    def test_fail(self):
        """ Test our ability to fail """

        result = ComputationalResult()
        result.save()

        sample = Sample()
        sample.accession_code = 'XXX'
        sample.title = 'XXX'
        sample.organism = Organism.get_object_for_name("HOMO_SAPIENS")
        sample.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        computed_file = ComputedFile()
        computed_file.filename = "NOT_REAL.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {'GSE51081': ['XXX']}
        ds.aggregate_by = 'EXPERIMENT'
        ds.scale_by = 'MINMAX'
        ds.email_address = "null@derp.com"
        ds.save()
        dsid = ds.id

        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)
        ds = Dataset.objects.get(id=dsid)
        print(ds.failure_reason)
        print(final_context['dataset'].failure_reason)
        self.assertFalse(ds.success)
        self.assertNotEqual(ds.failure_reason, "")

    @tag("smasher")
    def test_no_smash_all_diff_species(self):
        """ Smashing together with 'ALL' with different species should
        cause a 0 length data frame after inner join. """

        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51081"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        sample = Sample()
        sample.accession_code = 'GSM1237810'
        sample.title = 'GSM1237810'
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
        computed_file.filename = "GSM1237810_T09-1084.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        result = ComputationalResult()
        result.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51084"
        experiment.save()

        mus_mus = Organism.get_object_for_name("MUS_MUSCULUS")

        sample = Sample()
        sample.accession_code = 'GSM1238108'
        sample.title = 'GSM1238108'
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
        computed_file.filename = "GSM1238108-tbl-1.txt"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {'GSE51081': ['GSM1237810'], 'GSE51084': ['GSM1238108']}
        ds.aggregate_by = 'ALL'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        print("This is supposed to fail!:")
        print(ds.failure_reason)
        print(final_context['dataset'].failure_reason)

        self.assertFalse(ds.success)
        self.assertNotEqual(ds.failure_reason, "")
        self.assertEqual(len(final_context['merged']), 0)

        ds.delete()
        job.delete()

    @tag("smasher")
    def test_dualtech_smash(self):
        """ """

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        # MICROARRAY TECH
        experiment = Experiment()
        experiment.accession_code = "GSM1487313"
        experiment.save()

        result = ComputationalResult()
        result.save()

        gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS")

        sample = Sample()
        sample.accession_code = 'GSM1487313'
        sample.title = 'GSM1487313'
        sample.organism = gallus_gallus
        sample.technology="MICROARRAY"
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
        computed_file.filename = "GSM1487313_liver.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        # RNASEQ TECH
        experiment = Experiment()
        experiment.accession_code = "SRS332914"
        experiment.save()

        result2 = ComputationalResult()
        result2.save()

        sample = Sample()
        sample.accession_code = 'SRS332914'
        sample.title = 'SRS332914'
        sample.organism = gallus_gallus
        sample.technology="RNA-SEQ"
        sample.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result2
        sra.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = experiment
        esa.sample = sample
        esa.save()

        computed_file = ComputedFile()
        computed_file.filename = "SRP149598_gene_lengthScaledTPM.tsv"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result2
        computed_file.size_in_bytes = 234
        computed_file.is_smashable = True
        computed_file.save()

        # CROSS-SMASH BY SPECIES
        ds = Dataset()
        ds.data = {'GSM1487313': ['GSM1487313'], 'SRS332914': ['SRS332914']}
        ds.aggregate_by = 'SPECIES'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        self.assertTrue(ds.is_cross_technology())

        final_context = smasher.smash(pj.pk, upload=False)

        # THEN BY EXPERIMENT
        ds.aggregate_by = 'EXPERIMENT'
        ds.save()

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        final_context = smasher.smash(pj.pk, upload=False)

    @tag("smasher")
    def test_sanity_imports(self):
        """ Sci imports can be tricky, make sure this works. """

        import numpy
        import scipy
        import matplotlib
        import pandas
        import sklearn
        import sympy

import os
import shutil
import sys
import zipfile
from io import StringIO

from django.core.management import call_command
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
    SampleAnnotation,
    SampleResultAssociation,
    ExperimentSampleAssociation,
    Dataset,
    ProcessorJobDatasetAssociation,
    SampleComputedFileAssociation
)
from data_refinery_workers.processors import smasher


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

    sample_annotation = SampleAnnotation()
    sample_annotation.data = {'hi': 'friend'}
    sample_annotation.sample = sample
    sample_annotation.save()

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

    computed_file = ComputedFile()
    computed_file.filename = "GSM1237812_S97-PURE.DAT"
    computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
    computed_file.result = result
    computed_file.size_in_bytes = 123
    computed_file.is_smashable = False
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
    ds.quantile_normalize = False
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

        anno_samp = Sample.objects.get(accession_code='GSM1237810')
        self.assertTrue('hi' in anno_samp.to_metadata_dict()['annotations'][0].keys())

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

        # XXX: agg_type 'SPECIES' hangs on Linux, not OSX.
        # Don't know why yet.
        # for ag_type in ['ALL', 'EXPERIMENT', 'SPECIES']:
        for ag_type in ['ALL', 'EXPERIMENT']:
            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.aggregate_by = ag_type
            dataset.save()

            print ("Smashing " + ag_type)
            final_context = smasher.smash(job.pk, upload=False)
            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context['output_file']), 0)
            self.assertEqual(final_context['dataset'].is_processed, True)

            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.is_processed = False
            dataset.save()

            # Cleanup
            os.remove(final_context['output_file'])
            job.start_time = None
            job.end_time = None
            job.save()

        for scale_type in ['NONE', 'MINMAX', 'STANDARD', 'ROBUST']:
            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.aggregate_by = 'EXPERIMENT'
            dataset.scale_by = scale_type
            dataset.save()

            print ("Smashing " + scale_type)
            final_context = smasher.smash(job.pk, upload=False)
            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context['output_file']), 0)
            self.assertEqual(final_context['dataset'].is_processed, True)

            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.is_processed = False
            dataset.save()

            # Cleanup
            os.remove(final_context['output_file'])
            job.start_time = None
            job.end_time = None
            job.save()

        # Stats
        for scale_type in ['MINMAX', 'STANDARD']:
            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.aggregate_by = 'EXPERIMENT'
            dataset.scale_by = scale_type
            dataset.save()

            print("###")
            print("# " + scale_type)
            print('###')

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

            self.assertTrue('GSE51081_metadata.tsv' in namelist)
            self.assertTrue('metadata.json' in namelist)
            self.assertTrue('README.md' in namelist)
            self.assertTrue('LICENSE.TXT' in namelist)
            self.assertTrue('GSE51081.tsv' in namelist)

            os.remove(final_context['output_file'])
            job.start_time = None
            job.end_time = None
            job.save()


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
        ds.quantile_normalize = False
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
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)
        print(ds.failure_reason)
        print(final_context['dataset'].failure_reason)

        self.assertFalse(ds.success)
        self.assertNotEqual(ds.failure_reason, "")
        self.assertEqual(len(final_context['original_merged']), 0)

    @tag("smasher")
    def test_no_smash_dupe(self):
        """ """

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

        sample = Sample()
        sample.accession_code = 'GSM1237811'
        sample.title = 'GSM1237811'
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

        result = ComputationalResult()
        result.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {'GSE51081': ['GSM1237810', 'GSM1237811']}
        ds.aggregate_by = 'ALL'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        self.assertTrue(ds.success)
        for column in final_context['original_merged'].columns:
            self.assertTrue('_x' not in column)

    @tag("smasher")
    def test_log2(self):
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        # Has non-log2 data:
        # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44421
        # ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44421/miniml/GSE44421_family.xml.tgz
        experiment = Experiment()
        experiment.accession_code = "GSE44421"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        sample = Sample()
        sample.accession_code = 'GSM1084806'
        sample.title = 'GSM1084806'
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
        computed_file.filename = "GSM1084806-tbl-1.txt"
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
        sample.accession_code = 'GSM1084807'
        sample.title = 'GSM1084807'
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
        computed_file.filename = "GSM1084807-tbl-1.txt"
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
        ds.data = {'GSE44421': ['GSM1084806', 'GSM1084807']}
        ds.aggregate_by = 'EXPERIMENT'
        ds.scale_by = 'MINMAX'
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)

        self.assertTrue(final_context['success'])

    @tag("smasher")
    def test_dualtech_smash(self):
        """ """

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        # MICROARRAY TECH
        experiment = Experiment()
        experiment.accession_code = "GSE1487313"
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

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        # RNASEQ TECH
        experiment2 = Experiment()
        experiment2.accession_code = "SRS332914"
        experiment2.save()

        result2 = ComputationalResult()
        result2.save()

        sample2 = Sample()
        sample2.accession_code = 'SRS332914'
        sample2.title = 'SRS332914'
        sample2.organism = gallus_gallus
        sample2.technology = "RNA-SEQ"
        sample2.save()

        sra2 = SampleResultAssociation()
        sra2.sample = sample2
        sra2.result = result2
        sra2.save()

        esa2 = ExperimentSampleAssociation()
        esa2.experiment = experiment2
        esa2.sample = sample2
        esa2.save()

        computed_file2 = ComputedFile()
        computed_file2.filename = "SRP149598_gene_lengthScaledTPM.tsv"
        computed_file2.absolute_file_path = "/home/user/data_store/PCL/" + computed_file2.filename
        computed_file2.result = result2
        computed_file2.size_in_bytes = 234
        computed_file2.is_smashable = True
        computed_file2.save()

        assoc2 = SampleComputedFileAssociation()
        assoc2.sample = sample2
        assoc2.computed_file = computed_file2
        assoc2.save()

        # CROSS-SMASH BY SPECIES
        ds = Dataset()
        ds.data = {'GSE1487313': ['GSM1487313'], 'SRX332914': ['SRS332914']}
        ds.aggregate_by = 'SPECIES'
        ds.scale_by = 'STANDARD'
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        self.assertTrue(ds.is_cross_technology())
        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context['output_file']))
        os.remove(final_context['output_file'])
        self.assertEqual(len(final_context['final_frame'].columns), 2)

        # THEN BY EXPERIMENT
        ds.aggregate_by = 'EXPERIMENT'
        ds.save()

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        pj.start_time = None
        pj.end_time = None
        pj.save()

        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context['output_file']))
        os.remove(final_context['output_file'])
        self.assertEqual(len(final_context['final_frame'].columns), 1)

        # THEN BY ALL
        ds.aggregate_by = 'ALL'
        ds.save()

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        pj.start_time = None
        pj.end_time = None
        pj.save()
        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context['output_file']))
        self.assertEqual(len(final_context['final_frame'].columns), 2)

    @tag("smasher")
    def test_sanity_imports(self):
        """ Sci imports can be tricky, make sure this works. """

        import numpy
        import scipy
        import matplotlib
        import pandas
        import sklearn
        import sympy

    @tag("smasher")
    def test_get_synced_files(self):
        """ """
        result = ComputationalResult()
        result.save()

        computed_file = ComputedFile()
        computed_file.s3_key = "all_the_things.jpg"
        computed_file.s3_bucket = "data-refinery-test-assets"
        computed_file.filename = "all_the_things.jpg"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 9001
        computed_file.is_smashable = False
        computed_file.save()

        # Make sure it's not there
        try:
            os.remove("/home/user/data_store/PCL/" + computed_file.filename)
        except OSError:
            pass

        # We do this twice, once to get from S3 and once to get from local disk.
        afp = computed_file.get_synced_file_path(force=True)
        self.assertTrue(os.path.exists(afp))

        afp = computed_file.get_synced_file_path(force=True)
        self.assertTrue(os.path.exists(afp))

    @tag("smasher")
    def test_notify(self):

        ds = Dataset()
        ds.data = {'GSM1487313': ['GSM1487313'], 'SRS332914': ['SRS332914']}
        ds.aggregate_by = 'SPECIES'
        ds.scale_by = 'STANDARD'
        ds.email_address = "shoopdawoop@mailinator.com"
        ds.quantile_normalize = False
        ds.save()

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        job_context = {}
        job_context['job'] = pj
        job_context['dataset'] = ds
        job_context['upload'] = True
        job_context['result_url'] = 'https://s3.amazonaws.com/data-refinery-test-assets/all_the_things.jpg'

        final_context = smasher._notify(job_context)
        self.assertTrue(final_context.get('success', True))


class CompendiaTestCase(TestCase):
    """ Testing management commands are hard. Since there is always an explicit sys.exit (which is really an Exception),
    we have to do weird stdio rerouting to capture the result. Really, these are just sanity tests."""

    @tag("smasher")
    def test_call_create(self):
        old_stderr = sys.stderr
        old_stdout = sys.stdout
        csio_err = StringIO()
        csio_out = StringIO()
        sys.stderr = csio_err
        sys.stdout = csio_out
        self.assertRaises(BaseException, call_command, 'create_compendia')
        sys.stderr = old_stderr
        sys.stdout = old_stdout

    @tag("smasher")
    def test_fetch_create(self):
        old_stderr = sys.stderr
        old_stdout = sys.stdout
        csio_err = StringIO()
        csio_out = StringIO()
        sys.stderr = csio_err
        sys.stdout = csio_out
        self.assertRaises(BaseException, call_command, 'fetch_compendia')
        sys.stderr = old_stderr
        sys.stdout = old_stdout

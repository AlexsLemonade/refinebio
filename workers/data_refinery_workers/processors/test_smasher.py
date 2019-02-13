# -*- coding: utf-8 -*- 

import csv
import json
import os
import shutil
import sys
import zipfile
from io import StringIO

import pandas as pd

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
    SampleComputedFileAssociation,
    ComputationalResultAnnotation
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
        self.assertTrue('hi' in anno_samp.to_metadata_dict()['refinebio_annotations'][0].keys())

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

            self.assertFalse(True in final_frame.index.str.contains('AFFX-'))
            self.assertTrue('GSE51081/metadata_GSE51081.tsv' in namelist)
            self.assertTrue('aggregated_metadata.json' in namelist)
            self.assertTrue('README.md' in namelist)
            self.assertTrue('LICENSE.TXT' in namelist)
            self.assertTrue('GSE51081/GSE51081.tsv' in namelist)

            os.remove(final_context['output_file'])
            job.start_time = None
            job.end_time = None
            job.save()

        for scale_type in ['MINMAX', 'STANDARD']:
            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.aggregate_by = 'SPECIES'
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

            self.assertTrue('HOMO_SAPIENS/metadata_HOMO_SAPIENS.tsv' in namelist)
            self.assertTrue('aggregated_metadata.json' in namelist)
            self.assertTrue('README.md' in namelist)
            self.assertTrue('LICENSE.TXT' in namelist)
            self.assertTrue('HOMO_SAPIENS/HOMO_SAPIENS.tsv' in namelist)

            os.remove(final_context['output_file'])
            job.start_time = None
            job.end_time = None
            job.save()

        for scale_type in ['MINMAX', 'STANDARD']:
            dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
            dataset.aggregate_by = 'ALL'
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

            self.assertTrue('ALL/metadata_ALL.tsv' in namelist)
            self.assertTrue('aggregated_metadata.json' in namelist)
            self.assertTrue('README.md' in namelist)
            self.assertTrue('LICENSE.TXT' in namelist)
            self.assertTrue('ALL/ALL.tsv' in namelist)

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
        self.assertNotEqual(final_context['unsmashable_files'], [])

    @tag("smasher")
    def test_no_smash_all_diff_species(self):
        """ Smashing together with 'ALL' with different species is a really weird behavior. 
        This test isn't really testing a normal case, just make sure that it's marking the
        unsmashable files.
        """

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

        self.assertEqual(final_context['unsmashable_files'], ['GSM1238108'])

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
    def test_no_smash_dupe_two(self):
        """ Tests the SRP051449 case, where the titles collide. Also uses a real QN target file."""

        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        experiment = Experiment()
        experiment.accession_code = "SRP051449"
        experiment.save()

        result = ComputationalResult()
        result.save()

        danio_rerio = Organism.get_object_for_name("DANIO_RERIO")

        sample = Sample()
        sample.accession_code = 'SRR1731761'
        sample.title = 'Danio rerio'
        sample.organism = danio_rerio
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
        computed_file.filename = "SRR1731761_output_gene_lengthScaledTPM.tsv"
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
        sample.accession_code = 'SRR1731762'
        sample.title = 'Danio rerio'
        sample.organism = danio_rerio
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
        computed_file.filename = "SRR1731762_output_gene_lengthScaledTPM.tsv"
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

        ds = Dataset()
        ds.data = {'SRP051449': ['SRR1731761', 'SRR1731762']}
        ds.aggregate_by = 'SPECIES'
        ds.scale_by = 'NONE'
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = True
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        cr = ComputationalResult()
        cr.save()

        computed_file = ComputedFile()
        computed_file.filename = "danio_target.tsv"
        computed_file.absolute_file_path = "/home/user/data_store/PCL/" + computed_file.filename
        computed_file.result = cr
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = False
        computed_file.save()

        cra = ComputationalResultAnnotation()
        cra.data = {'organism_id': danio_rerio.id, 'is_qn': True}
        cra.result = cr
        cra.save()

        final_context = smasher.smash(job.pk, upload=False)
        self.assertTrue(final_context['success'])

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
        computed_file.sha1 = "36cf21c08d461f74ddb0f2edb6257afee309c4a4"
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
    """Testing management commands are hard.  Since there is always an explicit
    sys.exit (which is really an Exception), we have to do weird stdio rerouting
    to capture the result. Really, these are just sanity tests.
    """

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


class AggregationTestCase(TestCase):
    """Test the tsv file generation."""
    def setUp(self):
        self.metadata = {
            'experiments': {
                "E-GEOD-44719": {
                    "accession_code": "E-GEOD-44719",
                    "sample_titles": [ "IFNa DC_LB016_IFNa", "undefined_sample" ]
                }
            },

            'samples': {
                "IFNa DC_LB016_IFNa": {  # Sample #1 is an ArrayExpress sample
                    "refinebio_title": "IFNa DC_LB016_IFNa",
                    "refinebio_accession_code": "E-GEOD-44719-GSM1089311",
                    "refinebio_source_database": "ARRAY_EXPRESS",
                    "refinebio_organism": "fake_species",
                    ############# Annotations will be de-composed. #############
                    "refinebio_annotations": [
                        # annotation #1
                        {
                            "detected_platform": "illuminaHumanv3",
                            "detection_percentage": 98.44078,
                            "mapped_percentage": 100.0
                        },
                        # annotation #2
                        {
                            "assay": { "name": "GSM1089311" },
                            # Special field that will be taken out as separate columns
                            "characteristic": [
                                { "category": "cell population",
                                  "value": "IFNa DC"
                                },
                                { "category": "dose",   # also available in "variable"
                                  "value": "1 mL"
                                },
                                { "category": "donor id",
                                  "value": "LB016"
                                }
                            ],
                            # Another special field in Array Express sample
                            "variable": [
                                { "name": "dose",  # also available in "characteristic"
                                  "value": "1 mL"
                                },
                                { "name": "stimulation",
                                  "value": "IFNa"
                                }
                            ],
                            # "source" field in Array Express sample annotation will be
                            # skipped in tsv file.
                            'source': {
                                'name': 'GSM1288968 1',
                                'comment': [
                                    { 'name': 'Sample_source_name',
                                      'value': 'pineal glands at CT18, after light exposure'
                                    },
                                    { 'name': 'Sample_title',
                                      'value': 'Pineal_Light_CT18'
                                    }
                                ]
                            },

                            # For single-key object whose key is "name",
                            # the key will be ignored in tsv file.
                            "extract": { "name": "GSM1089311 extract 1" }
                        }
                    ]  # end of annotations
                },  # end of sample #1

                "Bone.Marrow_OA_No_ST03": {  # Sample #2 is a GEO sample
                    "refinebio_title": "Bone.Marrow_OA_No_ST03",
                    "refinebio_accession_code": "GSM1361050",
                    "refinebio_source_database": "GEO",
                    "refinebio_organism": "homo_sapiens",

                    "refinebio_annotations": [
                        {
                            "channel_count": [ "1" ],

                            # Special field that will be taken out as separate columns
                            "characteristics_ch1": [
                                "tissue: Bone Marrow",
                                "disease: OA",
                                "serum: Low Serum"
                            ],

                            # For single-element array, the element will
                            # be saved directly in tsv file.
                            "contact_address": [ "Crown Street" ],
                            "contact_country": [ "United Kingdom" ],
                            "data_processing": [ "Data was processed and normalized" ],
                            "geo_accession": [ "GSM1361050" ],
                        }
                    ]  # end of annotations
                }  # end of sample #2

            }  # end of "samples"
        }

        self.smash_path = "/tmp/"

    @tag("smasher")
    def test_columns(self):
        columns = smasher._get_tsv_columns(self.metadata['samples'])
        self.assertEqual(len(columns), 21)
        self.assertEqual(columns[0], 'refinebio_accession_code')
        self.assertTrue('refinebio_accession_code' in columns)
        self.assertTrue('cell population' in columns)
        self.assertTrue('dose' in columns)
        self.assertTrue('stimulation' in columns)
        self.assertTrue('serum' in columns)

    @tag("smasher")
    def test_all_samples(self):
        """Check tsv file that includes all sample metadata."""

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='ALL')
        }
        smasher._write_tsv_json(job_context, self.metadata, self.smash_path)
        tsv_filename = self.smash_path + "ALL/metadata_ALL.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))

        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                if row['refinebio_accession_code'] == 'E-GEOD-44719-GSM1089311':
                    self.assertEqual(row['cell population'], 'IFNa DC') # ArrayExpress specific
                    self.assertEqual(row['dose'], '1 mL')               # ArrayExpress specific
                    self.assertFalse('source' in row)                   # ArrayExpress specific
                    self.assertEqual(row['detection_percentage'], '98.44078')
                    self.assertEqual(row["extract"], "GSM1089311 extract 1")
                elif row['refinebio_accession_code'] == 'GSM1361050':
                    self.assertEqual(row['tissue'], 'Bone Marrow')      # GEO specific
                    self.assertEqual(row['refinebio_organism'], 'homo_sapiens')
                    self.assertEqual(row["contact_address"], "Crown Street")

        self.assertEqual(row_num, 1)  # only two data rows in tsv file
        os.remove(tsv_filename)

    @tag("smasher")
    def test_experiment(self):
        """Check tsv file that is aggregated by experiment."""

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='EXPERIMENT')
        }
        smasher._write_tsv_json(job_context, self.metadata, self.smash_path)

        tsv_filename = self.smash_path + "E-GEOD-44719/metadata_E-GEOD-44719.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))

        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                self.assertEqual(row['refinebio_accession_code'], 'E-GEOD-44719-GSM1089311')
                self.assertEqual(row['cell population'], 'IFNa DC')  # ArrayExpress specific
                self.assertEqual(row['dose'], '1 mL')                # ArrayExpress specific
                self.assertEqual(row['detection_percentage'], '98.44078')

        self.assertEqual(row_num, 0) # only one data row in tsv file
        os.remove(tsv_filename)

    @tag("smasher")
    def test_unicode_writer(self):
        self.unicode_metadata = {
            'experiments': {
                "E-GEOD-😎": {
                    "accession_code": "E-GEOD-😎",
                    "sample_titles": [ "😎", "undefined_sample" ]
                }
            },
            'samples': {
                "😎": {  # Sample #1 is an ArrayExpress sample
                    "refinebio_😎": "😎",
                    "refinebio_accession_code": "eyy",
                    "refinebio_annotations": [
                        # annotation #1
                        {
                            "😎😎": "😎😎",
                            "detection_percentage": 98.44078,
                            "mapped_percentage": 100.0
                        }
                    ]  # end of annotations
                },  # end of sample #1
            }  # end of "samples"
        }
        self.smash_path = "/tmp/"

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='ALL'),
            'input_files': {
            }
        }
        final_context = smasher._write_tsv_json(job_context, self.unicode_metadata, self.smash_path)
        reso = final_context[0]
        with open(reso, encoding='utf-8') as tsv_file: 
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                print(str(row).encode('utf-8'))

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='EXPERIMENT'),
            'input_files': {
            }
        }
        final_context = smasher._write_tsv_json(job_context, self.unicode_metadata, self.smash_path)
        reso = final_context[0]
        with open(reso, encoding='utf-8') as tsv_file: 
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                print(str(row).encode('utf-8'))

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='SPECIES'),
            'input_files': {
                'homo_sapiens': [], # only the key matters in this test
                'fake_species': []  # only the key matters in this test
            }
        }
        final_context = smasher._write_tsv_json(job_context, self.unicode_metadata, self.smash_path)
        reso = final_context[0]
        with open(reso, encoding='utf-8') as tsv_file: 
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                print(str(row).encode('utf-8'))

    @tag("smasher")
    def test_species(self):
        """Check tsv file that is aggregated by species."""

        job_context = {
            'dataset': Dataset.objects.create(aggregate_by='SPECIES'),
            'input_files': {
                'homo_sapiens': [], # only the key matters in this test
                'fake_species': []  # only the key matters in this test
            }
        }
        # Generate two TSV files, one should include only "GSM1361050",
        # and the other should include only "E-GEOD-44719-GSM1089311".
        smasher._write_tsv_json(job_context, self.metadata, self.smash_path)

        # Test tsv file of "homo_sapiens"
        tsv_filename = self.smash_path + "homo_sapiens/metadata_homo_sapiens.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))
        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                self.assertEqual(row['refinebio_accession_code'], 'GSM1361050')
                self.assertEqual(row['tissue'], 'Bone Marrow')      # GEO specific
                self.assertEqual(row['refinebio_organism'], 'homo_sapiens')

        self.assertEqual(row_num, 0) # only one data row in tsv file
        os.remove(tsv_filename)

        # Test json file of "homo_sapiens"
        json_filename = self.smash_path + "homo_sapiens/metadata_homo_sapiens.json"
        self.assertTrue(os.path.isfile(json_filename))
        with open(json_filename) as json_fp:
            species_metadada = json.load(json_fp)
        self.assertEqual(species_metadada['species'], 'homo_sapiens')
        self.assertEqual(len(species_metadada['samples']), 1)
        self.assertEqual(species_metadada['samples'][0]['refinebio_accession_code'],
                         'GSM1361050')
        #os.remove(json_filename)

        # Test tsv file of "fake_species"
        tsv_filename = self.smash_path + "fake_species/metadata_fake_species.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))
        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row_num, row in enumerate(reader):
                self.assertEqual(row['refinebio_accession_code'], 'E-GEOD-44719-GSM1089311')
                self.assertEqual(row['cell population'], 'IFNa DC')  # ArrayExpress specific
                self.assertEqual(row['dose'], '1 mL')                # ArrayExpress specific
                self.assertEqual(row['detection_percentage'], '98.44078')

        self.assertEqual(row_num, 0) # only one data row in tsv file
        os.remove(tsv_filename)

        # Test json file of "fake_species"
        json_filename = self.smash_path + "fake_species/metadata_fake_species.json"
        self.assertTrue(os.path.isfile(json_filename))

        with open(json_filename) as json_fp:
            species_metadada = json.load(json_fp)
        self.assertEqual(species_metadada['species'], 'fake_species')
        self.assertEqual(len(species_metadada['samples']), 1)
        self.assertEqual(species_metadada['samples'][0]['refinebio_accession_code'],
                         'E-GEOD-44719-GSM1089311')
        os.remove(json_filename)

    @tag("smasher")
    def test_bad_overlap(self):

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
        computed_file.filename = "big.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/BADSMASH/" + computed_file.filename
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
        computed_file.filename = "small.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/BADSMASH/" + computed_file.filename
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
        ds.aggregate_by = 'ALL' # [ALL or SPECIES or EXPERIMENT]
        ds.scale_by = 'NONE' # [NONE or MINMAX or STANDARD or ROBUST]
        ds.email_address = "null@derp.com"
        #ds.email_address = "miserlou+heyo@gmail.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        # Now, make sure the bad can't zero this out.
        sample = Sample()
        sample.accession_code = 'GSM999'
        sample.title = 'GSM999'
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
        computed_file.filename = "bad.PCL"
        computed_file.absolute_file_path = "/home/user/data_store/BADSMASH/" + computed_file.filename
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {'GSE51081': ['GSM1237810', 'GSM1237812', 'GSM999']}
        ds.aggregate_by = 'ALL' # [ALL or SPECIES or EXPERIMENT]
        ds.scale_by = 'NONE' # [NONE or MINMAX or STANDARD or ROBUST]
        ds.email_address = "null@derp.com"
        #ds.email_address = "miserlou+heyo@gmail.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)

        self.assertEqual(len(final_context['final_frame']), 4)

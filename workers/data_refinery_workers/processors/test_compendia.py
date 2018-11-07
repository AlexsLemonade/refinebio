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
    ProcessorJobDatasetAssociation,
    SampleAnnotation,
    SampleResultAssociation
)
from data_refinery_workers.processors import create_compendia, smasher, utils


class CompendiaTestCase(TestCase):

    @tag('compendia')
    def test_sanity(self):
        print("Hey!")

    @tag('compendia')
    def test_create_compendia(self):
        job = ProcessorJob()
        job.pipeline_applied = "COMPENDIA"
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

        dset = Dataset()
        dset.data = {'GSE51081': ['GSM1237810', 'GSM1237812']}
        dset.scale_by = 'NONE'
        dset.aggregate_by = 'SPECIES'
        dset.quantile_normalize = False
        dset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dset
        pjda.save()

        final_context = create_compendia.create_compendia(job.id)

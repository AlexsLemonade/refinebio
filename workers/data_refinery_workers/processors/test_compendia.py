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

        dset = Dataset()
        dset.data = {'GSE1487313': ['GSM1487313'], 'SRX332914': ['SRS332914']}
        dset.scale_by = 'NONE'
        dset.aggregate_by = 'SPECIES'
        dset.quantile_normalize = False
        dset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dset
        pjda.save()

        final_context = create_compendia.create_compendia(job.id)

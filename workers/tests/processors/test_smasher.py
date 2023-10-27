# -*- coding: utf-8 -*-

import csv
import json
import os
import sys
import zipfile
from io import StringIO
from unittest.mock import MagicMock, patch

from django.core.management import call_command
from django.test import TransactionTestCase, tag

import pandas as pd
import vcr
from tests.utils import ProcessorJobTestCaseMixin

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Contribution,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    OntologyTerm,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SampleAnnotation,
    SampleAttribute,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_workers.processors import smasher, smashing_utils


def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SMASHER"
    pj.save()

    experiment = Experiment()
    experiment.accession_code = "GSE51081"
    experiment.save()

    result = ComputationalResult()
    result.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=9606)

    sample = Sample()
    sample.accession_code = "GSM1237810"
    sample.title = "GSM1237810"
    sample.organism = homo_sapiens
    sample.save()

    sample_annotation = SampleAnnotation()
    sample_annotation.data = {"hi": "friend"}
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
    computed_file.filename = "SRP149598_gene_lengthScaledTPM.tsv"
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
    sample.accession_code = "GSM1237812"
    sample.title = "GSM1237812"
    sample.organism = homo_sapiens
    sample.save()

    esa = ExperimentSampleAssociation()
    esa.experiment = experiment
    esa.sample = sample
    esa.save()

    result = ComputationalResult()
    result.save()

    sra = SampleResultAssociation()
    sra.sample = sample
    sra.result = result
    sra.save()

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
    ds.data = {"GSE51081": ["GSM1237810", "GSM1237812"]}
    ds.aggregate_by = "EXPERIMENT"
    ds.scale_by = "STANDARD"
    ds.email_address = "null@derp.com"
    ds.quantile_normalize = False
    ds.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = pj
    pjda.dataset = ds
    pjda.save()

    return pj


def prepare_dual_tech_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SMASHER"
    pj.save()

    # MICROARRAY TECH
    experiment = Experiment()
    experiment.accession_code = "GSE1487313"
    experiment.save()

    result = ComputationalResult()
    result.save()

    gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS", taxonomy_id=9031)

    sample = Sample()
    sample.accession_code = "GSM1487313"
    sample.title = "GSM1487313"
    sample.organism = gallus_gallus
    sample.technology = "MICROARRAY"
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
    experiment2.accession_code = "SRP332914"
    experiment2.save()

    result2 = ComputationalResult()
    result2.save()

    sample2 = Sample()
    sample2.accession_code = "SRR332914"
    sample2.title = "SRR332914"
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
    return pj


def _to_accession_code_dict(d):
    return {k: set(map(lambda x: x.accession_code, v)) for k, v in d.items()}


class SmasherTestCase(TransactionTestCase, ProcessorJobTestCaseMixin):
    def _test_all_scale_types(self, ag_type, ag_specific_tests):
        job = prepare_job()

        anno_samp = Sample.objects.get(accession_code="GSM1237810")
        self.assertTrue("hi" in anno_samp.to_metadata_dict()["refinebio_annotations"][0].keys())

        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)
        dataset = Dataset.objects.filter(id__in=relations.values("dataset_id")).first()
        job_context_check = {}
        job_context_check["dataset"] = dataset
        job_context_check["samples"] = dataset.get_samples()
        job_context_check["experiments"] = dataset.get_experiments()
        self.assertEqual(len(job_context_check["samples"]), 2)
        self.assertEqual(len(job_context_check["experiments"]), 1)

        print(f"\n######\n### {ag_type}\n######")

        for scale_type in ["NONE", "MINMAX", "STANDARD", "ROBUST"]:
            dataset = Dataset.objects.filter(id__in=relations.values("dataset_id")).first()
            dataset.aggregate_by = ag_type
            dataset.scale_by = scale_type
            dataset.save()

            print(f"\n###### Smashing with scale type {scale_type}\n")

            final_context = smasher.smash(job.pk, upload=False)
            # Make sure the file exists and is a valid size
            self.assertNotEqual(os.path.getsize(final_context["output_file"]), 0)
            self.assertEqual(final_context["dataset"].is_processed, True)

            final_frame = final_context.get("final_frame")

            # Sanity test that these frames can be computed upon
            final_frame.mean(axis=1)
            final_frame.min(axis=1)
            final_frame.max(axis=1)
            final_frame.std(axis=1)
            final_frame.median(axis=1)

            zf = zipfile.ZipFile(final_context["output_file"])

            ag_specific_tests({"final_frame": final_frame, "zipfile": zf})

            dataset = Dataset.objects.filter(id__in=relations.values("dataset_id")).first()
            dataset.is_processed = False
            dataset.save()

            os.remove(final_context["output_file"])
            job.start_time = None
            job.end_time = None
            job.save()

    @tag("smasher")
    def test_smasher_experiment_aggregation(self):
        def experiment_tests(args):
            final_frame = args["final_frame"]
            zf = args["zipfile"]
            namelist = zf.namelist()

            self.assertNotIn(True, final_frame.index.str.contains("AFFX-"))
            self.assertIn("GSE51081/metadata_GSE51081.tsv", namelist)
            self.assertIn("aggregated_metadata.json", namelist)
            self.assertIn("GSE51081/GSE51081.tsv", namelist)
            with zf.open("GSE51081/GSE51081.tsv") as f:
                tsv = pd.read_csv(f, delimiter="\t")
                # Check that the columns match, and that we filtered one of the files
                self.assertEqual(list(tsv), ["Gene", "GSM1237810", "GSM1237812"])

            self.assertIn("README.md", namelist)
            with open("/home/user/README_DATASET.md", "r") as ex, zf.open("README.md") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

            self.assertIn("LICENSE.TXT", namelist)
            with open("/home/user/LICENSE_DATASET.txt", "r") as ex, zf.open("LICENSE.TXT") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

        self._test_all_scale_types("EXPERIMENT", experiment_tests)

    @tag("smasher")
    def test_get_samples_by_experiment(self):
        job = prepare_job()
        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)
        dataset = Dataset.objects.filter(id__in=relations.values("dataset_id")).first()

        self.assertEqual(
            _to_accession_code_dict(dataset.get_samples_by_experiment()),
            {"GSE51081": {"GSM1237810", "GSM1237812"}},
        )

        dataset.aggregate_by = "EXPERIMENT"
        self.assertEqual(
            _to_accession_code_dict(dataset.get_aggregated_samples()),
            _to_accession_code_dict(dataset.get_samples_by_experiment()),
        )

        job = prepare_dual_tech_job()
        dataset = Dataset()
        dataset.data = {"GSE1487313": ["GSM1487313"], "SRP332914": ["SRR332914"]}
        dataset.aggregate_by = "EXPERIMENT"
        dataset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        self.assertEqual(
            _to_accession_code_dict(dataset.get_samples_by_experiment()),
            {"GSE1487313": {"GSM1487313"}, "SRP332914": {"SRR332914"}},
        )

        dataset.aggregate_by = "EXPERIMENT"
        self.assertEqual(
            _to_accession_code_dict(dataset.get_aggregated_samples()),
            _to_accession_code_dict(dataset.get_samples_by_experiment()),
        )

    @tag("smasher")
    def test_smasher_species_aggregation(self):
        def species_tests(args):
            final_frame = args["final_frame"]
            zf = args["zipfile"]
            namelist = zf.namelist()

            self.assertIn("HOMO_SAPIENS/metadata_HOMO_SAPIENS.tsv", namelist)
            self.assertIn("aggregated_metadata.json", namelist)
            self.assertIn("HOMO_SAPIENS/HOMO_SAPIENS.tsv", namelist)
            with zf.open("HOMO_SAPIENS/HOMO_SAPIENS.tsv") as f:
                tsv = pd.read_csv(f, delimiter="\t")
                # Check that the columns match, and that we filtered one of the files
                self.assertEqual(list(tsv), ["Gene", "GSM1237810", "GSM1237812"])

            self.assertIn("README.md", namelist)
            with open("/home/user/README_DATASET.md", "r") as ex, zf.open("README.md") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

            self.assertIn("LICENSE.TXT", namelist)
            with open("/home/user/LICENSE_DATASET.txt", "r") as ex, zf.open("LICENSE.TXT") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

        self._test_all_scale_types("SPECIES", species_tests)

    @tag("smasher")
    def test_get_samples_by_species(self):
        job = prepare_job()
        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)
        dataset = Dataset.objects.filter(id__in=relations.values("dataset_id")).first()
        dataset.aggregate_by = "SPECIES"
        dataset.save()

        self.assertEqual(
            _to_accession_code_dict(dataset.get_samples_by_species()),
            {"HOMO_SAPIENS": {"GSM1237810", "GSM1237812"}},
        )

        self.assertEqual(dataset.aggregate_by, "SPECIES")
        self.assertEqual(
            _to_accession_code_dict(dataset.get_aggregated_samples()),
            _to_accession_code_dict(dataset.get_samples_by_species()),
        )

        job = prepare_dual_tech_job()
        dataset = Dataset()
        dataset.data = {"GSE1487313": ["GSM1487313"], "SRP332914": ["SRR332914"]}
        dataset.aggregate_by = "SPECIES"
        dataset.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = dataset
        pjda.save()

        self.assertEqual(
            _to_accession_code_dict(dataset.get_samples_by_species()),
            {"GALLUS_GALLUS": {"GSM1487313", "SRR332914"}},
        )

        dataset.aggregate_by = "SPECIES"
        self.assertEqual(
            _to_accession_code_dict(dataset.get_aggregated_samples()),
            _to_accession_code_dict(dataset.get_samples_by_species()),
        )

    @tag("smasher")
    def test_smasher_all_aggregation(self):
        def all_tests(args):
            final_frame = args["final_frame"]
            zf = args["zipfile"]
            namelist = zf.namelist()

            self.assertIn("ALL/metadata_ALL.tsv", namelist)
            self.assertIn("aggregated_metadata.json", namelist)
            self.assertIn("ALL/ALL.tsv", namelist)
            with zf.open("ALL/ALL.tsv") as f:
                tsv = pd.read_csv(f, delimiter="\t")
                # Check that the columns match, and that we filtered one of the files
                self.assertEqual(list(tsv), ["Gene", "GSM1237810", "GSM1237812"])

            self.assertIn("README.md", namelist)
            with open("/home/user/README_DATASET.md", "r") as ex, zf.open("README.md") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

            self.assertIn("LICENSE.TXT", namelist)
            with open("/home/user/LICENSE_DATASET.txt", "r") as ex, zf.open("LICENSE.TXT") as ac:
                self.assertEqual(ex.read(), ac.read().decode())

            with zf.open("aggregated_metadata.json") as aggregated_metadata:
                metadata_str = aggregated_metadata.read().decode()
                parsed_metadata = json.loads(metadata_str)
                # This dataset isn't quantile normalized, but we
                # should still be providing this value as False
                self.assertFalse(parsed_metadata["quantile_normalized"])

        self._test_all_scale_types("ALL", all_tests)

    @tag("smasher")
    def test_get_results(self):
        """Test our ability to collect the appropriate samples."""

        sample = Sample()
        sample.accession_code = "GSM45588"
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
        """Test our ability to fail"""

        result = ComputationalResult()
        result.save()

        sample = Sample()
        sample.accession_code = "XXX"
        sample.title = "XXX"
        sample.organism = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)
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
        ds.data = {"GSE51081": ["XXX"]}
        ds.aggregate_by = "EXPERIMENT"
        ds.scale_by = "MINMAX"
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
        print(final_context["dataset"].failure_reason)
        self.assertNotEqual(final_context["unsmashable_files"], [])

    @tag("smasher")
    def test_no_smash_all_diff_species(self):
        """Smashing together with 'ALL' with different species is a really weird behavior.
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

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "GSM1237810"
        sample.title = "GSM1237810"
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

        sample = Sample()
        sample.accession_code = "GSM1238108"
        sample.title = "GSM1238108"
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
        ds.data = {"GSE51081": ["GSM1237810"], "GSE51084": ["GSM1238108"]}
        ds.aggregate_by = "ALL"
        ds.scale_by = "STANDARD"
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
        print(final_context["dataset"].failure_reason)

        self.assertEqual(final_context["unsmashable_files"], ["GSM1237810_T09-1084.PCL"])

    @tag("smasher")
    def test_qn_targets_only(self):
        """ """
        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51088"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "GSM1237818"
        sample.title = "GSM1237818"
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
        computed_file.sha1 = (
            "08c7ea90b66b52f7cd9d9a569717a1f5f3874967"  # this matches with the downloaded file
        )
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
        ds.data = {"GSE51088": ["GSM1237818"]}
        ds.aggregate_by = "SPECIES"
        ds.scale_by = "STANDARD"
        ds.email_address = "null@derp.com"
        ds.quant_sf_only = True  # Make the dataset include quant.sf files only
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)

        # Check that the sample was really generated
        self.assertTrue(
            os.path.exists(final_context["output_dir"] + "/HOMO_SAPIENS/GSM1237818_quant.sf")
        )
        self.assertTrue(final_context["metadata"]["quant_sf_only"])
        self.assertEqual(final_context["metadata"]["num_samples"], 1)
        self.assertEqual(final_context["metadata"]["num_experiments"], 1)
        self.assertTrue("aggregate_by" not in final_context["metadata"])

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

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "GSM1237810"
        sample.title = "GSM1237810"
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
        sample.accession_code = "GSM1237811"
        sample.title = "GSM1237811"
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
        ds.data = {"GSE51081": ["GSM1237810", "GSM1237811"]}
        ds.aggregate_by = "ALL"
        ds.scale_by = "STANDARD"
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
        for column in final_context["original_merged"].columns:
            self.assertTrue("_x" not in column)

    @tag("smasher")
    def test_no_smash_dupe_two(self):
        """Tests the SRP051449 case, where the titles collide. Also uses a real QN target file."""

        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        experiment = Experiment()
        experiment.accession_code = "SRP051449"
        experiment.save()

        result = ComputationalResult()
        result.save()

        danio_rerio = Organism.get_object_for_name("DANIO_RERIO", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "SRR1731761"
        sample.title = "Danio rerio"
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
        sample.accession_code = "SRR1731762"
        sample.title = "Danio rerio"
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
        ds.data = {"SRP051449": ["SRR1731761", "SRR1731762"]}
        ds.aggregate_by = "SPECIES"
        ds.scale_by = "NONE"
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
        cra.data = {"organism_id": danio_rerio.id, "is_qn": True}
        cra.result = cr
        cra.save()

        danio_rerio.qn_target = cr
        danio_rerio.save()

        final_context = smasher.smash(job.pk, upload=False)
        self.assertSucceeded(job)

        # Test single file smash

        job = ProcessorJob()
        job.pipeline_applied = "SMASHER"
        job.save()

        ds = Dataset()
        ds.data = {"SRP051449": ["SRR1731761"]}
        ds.aggregate_by = "EXPERIMENT"
        ds.scale_by = "NONE"
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = True
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = job
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(job.pk, upload=False)
        self.assertSucceeded(job)

        zf = zipfile.ZipFile(final_context["output_file"])

        with zf.open("aggregated_metadata.json") as aggregated_metadata:
            metadata_str = aggregated_metadata.read().decode()
            parsed_metadata = json.loads(metadata_str)
            self.assertTrue(parsed_metadata["quantile_normalized"])

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

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)

        sample = Sample()
        sample.accession_code = "GSM1084806"
        sample.title = "GSM1084806"
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
        sample.accession_code = "GSM1084807"
        sample.title = "GSM1084807"
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
        ds.data = {"GSE44421": ["GSM1084806", "GSM1084807"]}
        ds.aggregate_by = "EXPERIMENT"
        ds.scale_by = "MINMAX"
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)

        self.assertSucceeded(pj)

    @tag("smasher")
    def test_dualtech_smash(self):
        """ """

        pj = prepare_dual_tech_job()
        # CROSS-SMASH BY SPECIES
        ds = Dataset()
        ds.data = {"GSE1487313": ["GSM1487313"], "SRP332914": ["SRR332914"]}
        ds.aggregate_by = "SPECIES"
        ds.scale_by = "STANDARD"
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        self.assertTrue(ds.is_cross_technology())
        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context["output_file"]))
        os.remove(final_context["output_file"])
        self.assertEqual(len(final_context["final_frame"].columns), 2)

        # THEN BY EXPERIMENT
        ds.aggregate_by = "EXPERIMENT"
        ds.save()

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        pj.start_time = None
        pj.end_time = None
        pj.save()

        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context["output_file"]))
        os.remove(final_context["output_file"])
        self.assertEqual(len(final_context["final_frame"].columns), 1)

        # THEN BY ALL
        ds.aggregate_by = "ALL"
        ds.save()

        dsid = ds.id
        ds = Dataset.objects.get(id=dsid)

        pj.start_time = None
        pj.end_time = None
        pj.save()
        final_context = smasher.smash(pj.pk, upload=False)
        self.assertTrue(os.path.exists(final_context["output_file"]))
        self.assertEqual(len(final_context["final_frame"].columns), 2)

    @tag("smasher")
    def test_quant_sf_only_smash(self):
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51088"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=9606)

        sample = Sample()
        sample.accession_code = "GSM1237818"
        sample.title = "GSM1237818"
        sample.organism = homo_sapiens
        sample.technology = "RNA-SEQ"
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
        computed_file.is_smashable = (
            False  # We should be able to handle quant.sf files that are not smashable (#2252)
        )
        computed_file.size_in_bytes = 123123
        computed_file.sha1 = (
            "08c7ea90b66b52f7cd9d9a569717a1f5f3874967"  # this matches with the downloaded file
        )
        computed_file.save()

        # Create a second sample whose quant.sf file is a real truncated file
        # from prod, and make sure that we filter it out properly
        result = ComputationalResult()
        result.save()

        sample = Sample()
        sample.accession_code = "GSM1237819"
        sample.title = "GSM1237819"
        sample.organism = homo_sapiens
        sample.technology = "RNA-SEQ"
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
        computed_file.s3_key = "smasher-test-truncated-quant.sf"
        computed_file.s3_bucket = "data-refinery-test-assets"
        computed_file.filename = "quant.sf"
        computed_file.absolute_file_path = (
            "/home/user/data_store/QUANT/smasher-test-truncated-quant.sf"
        )
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.size_in_bytes = 123123
        computed_file.sha1 = (
            "7a610039885be5f56b8ba29cf58d6555b1707ca5"  # this matches with the downloaded file
        )
        computed_file.save()

        ds = Dataset()
        ds.data = {"GSE51088": ["GSM1237818", "GSM1237819"]}
        ds.aggregate_by = "EXPERIMENT"
        ds.scale_by = "STANDARD"
        ds.email_address = "null@derp.com"
        ds.quant_sf_only = True  # Make the dataset include quant.sf files only
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)

        self.assertTrue(
            os.path.exists(final_context["output_dir"] + "/GSE51088/GSM1237818_quant.sf")
        )

        self.assertTrue(final_context["metadata"]["quant_sf_only"])
        self.assertEqual(final_context["metadata"]["num_samples"], 1)
        self.assertEqual(final_context["metadata"]["num_experiments"], 1)
        self.assertEqual(
            final_context["filtered_samples"]["GSM1237819"]["reason"],
            "This sample's quant.sf file was truncated and missing its header",
        )

    @tag("smasher")
    def test_sanity_imports(self):
        """Sci imports can be tricky, make sure this works."""

        import matplotlib  # noqa
        import numpy  # noqa
        import pandas  # noqa
        import scipy  # noqa
        import sklearn  # noqa
        import sympy  # noqa

    @tag("smasher")
    @vcr.use_cassette("/home/user/data_store/cassettes/smasher.get_synced_files.yaml")
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
    @patch("boto3.client")
    # We need to patch RUNNING_IN_CLOUD to get the smasher to attempt to send an email
    @patch("data_refinery_workers.processors.smasher.settings.RUNNING_IN_CLOUD", new=True)
    @patch("data_refinery_workers.processors.smasher.AWS_REGION", new="dummy region")
    def test_notify(self, mock_boto_client):
        ses_client_mock = MagicMock()
        mock_boto_client.return_value = ses_client_mock

        ds = Dataset()
        ds.data = {"GSM1487313": ["GSM1487313"], "SRS332914": ["SRS332914"]}
        ds.aggregate_by = "SPECIES"
        ds.scale_by = "STANDARD"
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
        job_context["job"] = pj
        job_context["dataset"] = ds
        job_context["upload"] = True
        job_context[
            "result_url"
        ] = "https://s3.amazonaws.com/data-refinery-test-assets/all_the_things.jpg"

        final_context = smasher._notify(job_context)
        self.assertTrue(final_context.get("success", True))
        mock_boto_client.assert_called_once_with("ses", region_name="dummy region")
        ses_client_mock.send_email.assert_called_once()

        # Now try explicitly setting notify_me true
        ses_client_mock.reset_mock()
        mock_boto_client.reset_mock()

        ds.notify_me = True
        ds.save()

        final_context = smasher._notify(job_context)

        self.assertTrue(final_context.get("success", True))
        mock_boto_client.assert_called_once_with("ses", region_name="dummy region")
        ses_client_mock.send_email.assert_called_once()
        self.assertEqual(
            ses_client_mock.send_email.call_args.kwargs["Source"],
            "refine.bio Datasets <noreply@refine.bio>",
        )

        # Now try explicitly setting notify_me false
        mock_boto_client.reset_mock()
        ses_client_mock.send_email.reset_mock()

        ds.notify_me = False
        ds.save()

        final_context = smasher._notify(job_context)

        self.assertTrue(final_context.get("success", True))
        mock_boto_client.assert_not_called()
        ses_client_mock.send_email.assert_not_called()


class CompendiaTestCase(TransactionTestCase):
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
        self.assertRaises(BaseException, call_command, "create_compendia")
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
        self.assertRaises(BaseException, call_command, "fetch_compendia")
        sys.stderr = old_stderr
        sys.stdout = old_stdout


class AggregationTestCase(TransactionTestCase):
    """Test the tsv file generation."""

    def setUp(self):
        self.metadata = {
            "experiments": {
                "E-GEOD-44719": {
                    "accession_code": "E-GEOD-44719",
                    "sample_accession_codes": [
                        "IFNa DC_LB016_IFNa",
                        "undefined_sample",
                    ],
                }
            },
            "samples": {
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
                            "mapped_percentage": 100.0,
                        },
                        # annotation #2
                        {
                            "assay": {"name": "GSM1089311"},
                            # Special field that will be taken out as separate columns
                            "characteristic": [
                                {"category": "cell population", "value": "IFNa DC"},
                                {
                                    "category": "dose",  # also available in "variable"
                                    "value": "1 mL",
                                },
                                {"category": "donor id", "value": "LB016"},
                            ],
                            # Another special field in Array Express sample
                            "variable": [
                                {
                                    "name": "dose",  # also available in "characteristic"
                                    "value": "1 mL",
                                },
                                {"name": "stimulation", "value": "IFNa"},
                            ],
                            # "source" field in Array Express sample annotation will be
                            # skipped in tsv file.
                            "source": {
                                "name": "GSM1288968 1",
                                "comment": [
                                    {
                                        "name": "Sample_source_name",
                                        "value": "pineal glands at CT18, after light exposure",
                                    },
                                    {
                                        "name": "Sample_title",
                                        "value": "Pineal_Light_CT18",
                                    },
                                ],
                            },
                            # For single-key object whose key is "name",
                            # the key will be ignored in tsv file.
                            "extract": {"name": "GSM1089311 extract 1"},
                        },
                    ],  # end of annotations
                },  # end of sample #1
                "Bone.Marrow_OA_No_ST03": {  # Sample #2 is a GEO sample
                    "refinebio_title": "Bone.Marrow_OA_No_ST03",
                    "refinebio_accession_code": "GSM1361050",
                    "refinebio_source_database": "GEO",
                    "refinebio_organism": "homo_sapiens",
                    "refinebio_annotations": [
                        {
                            "channel_count": ["1"],
                            # Special field that will be taken out as separate columns
                            "characteristics_ch1": [
                                "tissue: Bone Marrow",
                                "disease: OA",
                                "serum: Low Serum",
                            ],
                            # For single-element array, the element will
                            # be saved directly in tsv file.
                            "contact_address": ["Crown Street"],
                            "contact_country": ["United Kingdom"],
                            "data_processing": ["Data was processed and normalized"],
                            "geo_accession": ["GSM1361050"],
                        }
                    ],  # end of annotations
                },  # end of sample #2
            },  # end of "samples"
        }

        # Create a sample attribute to check for
        contribution, _ = Contribution.objects.get_or_create(
            source_name="MetaSRA",
            methods_url="https://pubmed.ncbi.nlm.nih.gov/28535296/",
        )

        s1 = Sample()
        s1.accession_code = "E-GEOD-44719-GSM1089311"
        s1.save()

        s2 = Sample()
        s2.accession_code = "GSM1361050"
        s2.save()

        age, _ = OntologyTerm.objects.get_or_create(
            ontology_term="EFO:0000246", human_readable_name="age"
        )

        attrib, _ = SampleAttribute.objects.get_or_create(name=age, source=contribution, sample=s1)
        attrib.set_value(3.0)
        attrib.save()

        self.smash_path = "/tmp/"

    @tag("smasher")
    def test_columns(self):
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        job_context = {"job": pj}

        columns = smashing_utils.get_tsv_columns(self.metadata["samples"])
        self.assertEqual(len(columns), 23)
        self.assertEqual(columns[0], "refinebio_accession_code")
        self.assertTrue("refinebio_accession_code" in columns)
        self.assertTrue("characteristic_cell population" in columns)
        self.assertTrue("characteristic_dose" in columns)
        self.assertTrue("characteristic_stimulation" in columns)
        self.assertTrue("characteristics_ch1_serum" in columns)
        self.assertTrue("MetaSRA_age" in columns)

    @tag("smasher")
    def test_all_samples(self):
        """Check tsv file that includes all sample metadata."""
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.metadata,
            "dataset": Dataset.objects.create(
                aggregate_by="ALL",
                data={
                    "GSE56409": ["GSM1361050"],
                    "E-GEOD-44719": ["E-GEOD-44719-GSM1089311"],
                },
            ),
        }
        smashing_utils.write_tsv_json(job_context)
        tsv_filename = self.smash_path + "ALL/metadata_ALL.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))

        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                if row["refinebio_accession_code"] == "E-GEOD-44719-GSM1089311":
                    self.assertEqual(
                        row["characteristic_cell population"], "IFNa DC"
                    )  # ArrayExpress specific
                    self.assertEqual(row["characteristic_dose"], "1 mL")  # ArrayExpress specific
                    self.assertFalse("source" in row)  # ArrayExpress specific
                    self.assertEqual(row["detection_percentage"], "98.44078")
                    self.assertEqual(row["extract"], "GSM1089311 extract 1")
                    self.assertEqual(row["experiment_accession"], "E-GEOD-44719")
                    self.assertEqual(row["MetaSRA_age"], "3.0")
                elif row["refinebio_accession_code"] == "GSM1361050":
                    self.assertEqual(
                        row["characteristics_ch1_tissue"], "Bone Marrow"
                    )  # GEO specific
                    self.assertEqual(row["refinebio_organism"], "homo_sapiens")
                    self.assertEqual(row["contact_address"], "Crown Street")
                    self.assertEqual(row["experiment_accession"], "GSE56409")

        self.assertEqual(row_num, 1)  # only two data rows in tsv file
        os.remove(tsv_filename)

    @tag("smasher")
    def test_experiment(self):
        """Check tsv file that is aggregated by experiment."""
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.metadata,
            "dataset": Dataset.objects.create(
                aggregate_by="EXPERIMENT",
                data={
                    "GSE56409": ["GSM1361050"],
                    "E-GEOD-44719": ["E-GEOD-44719-GSM1089311"],
                },
            ),
        }
        smashing_utils.write_tsv_json(job_context)

        tsv_filename = self.smash_path + "E-GEOD-44719/metadata_E-GEOD-44719.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))

        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                self.assertEqual(row["refinebio_accession_code"], "E-GEOD-44719-GSM1089311")
                self.assertEqual(row["experiment_accession"], "E-GEOD-44719")
                self.assertEqual(
                    row["characteristic_cell population"], "IFNa DC"
                )  # ArrayExpress specific
                self.assertEqual(row["characteristic_dose"], "1 mL")  # ArrayExpress specific
                self.assertEqual(row["detection_percentage"], "98.44078")
                self.assertEqual(row["MetaSRA_age"], "3.0")

        self.assertEqual(row_num, 0)  # only one data row in tsv file
        os.remove(tsv_filename)

    @tag("smasher")
    def test_unicode_writer(self):
        self.unicode_metadata = {
            "experiments": {
                "E-GEOD-😎": {
                    "accession_code": "E-GEOD-😎",
                    "sample_accession_codes": ["😎", "undefined_sample"],
                }
            },
            "samples": {
                "😎": {  # Sample #1 is an ArrayExpress sample
                    "refinebio_😎": "😎",
                    "refinebio_accession_code": "eyy",
                    "refinebio_annotations": [
                        # annotation #1
                        {
                            "😎😎": "😎😎",
                            "detection_percentage": 98.44078,
                            "mapped_percentage": 100.0,
                        }
                    ],  # end of annotations
                },  # end of sample #1
            },  # end of "samples"
        }
        self.smash_path = "/tmp/"

        s = Sample.objects.create(accession_code="eyy")

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.unicode_metadata,
            "dataset": Dataset.objects.create(aggregate_by="ALL"),
            "group_by_keys": [],
            "input_files": {},
        }
        final_context = smashing_utils.write_tsv_json(job_context)
        reso = final_context[0]
        with open(reso, encoding="utf-8") as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                print(str(row).encode("utf-8"))

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.unicode_metadata,
            "group_by_keys": [],
            "dataset": Dataset.objects.create(aggregate_by="EXPERIMENT"),
            "input_files": {},
        }
        final_context = smashing_utils.write_tsv_json(job_context)
        reso = final_context[0]
        with open(reso, encoding="utf-8") as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                print(str(row).encode("utf-8"))

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.unicode_metadata,
            "dataset": Dataset.objects.create(aggregate_by="SPECIES"),
            "group_by_keys": ["homo_sapiens", "fake_species"],
            "input_files": {
                "homo_sapiens": [],  # only the key matters in this test
                "fake_species": [],  # only the key matters in this test
            },
        }
        final_context = smashing_utils.write_tsv_json(job_context)
        reso = final_context[0]
        with open(reso, encoding="utf-8") as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                print(str(row).encode("utf-8"))

    @tag("smasher")
    def test_species(self):
        """Check tsv file that is aggregated by species."""

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        job_context = {
            "job": pj,
            "output_dir": self.smash_path,
            "metadata": self.metadata,
            "dataset": Dataset.objects.create(
                aggregate_by="SPECIES",
                data={
                    "GSE56409": ["GSM1361050"],
                    "E-GEOD-44719": ["E-GEOD-44719-GSM1089311"],
                },
            ),
            "group_by_keys": ["homo_sapiens", "fake_species"],
            "input_files": {
                "homo_sapiens": [],  # only the key matters in this test
                "fake_species": [],  # only the key matters in this test
            },
        }
        # Generate two TSV files, one should include only "GSM1361050",
        # and the other should include only "E-GEOD-44719-GSM1089311".
        smashing_utils.write_tsv_json(job_context)

        # Test tsv file of "homo_sapiens"
        tsv_filename = self.smash_path + "homo_sapiens/metadata_homo_sapiens.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))
        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                self.assertEqual(row["refinebio_accession_code"], "GSM1361050")
                self.assertEqual(row["experiment_accession"], "GSE56409")
                self.assertEqual(row["characteristics_ch1_tissue"], "Bone Marrow")  # GEO specific
                self.assertEqual(row["refinebio_organism"], "homo_sapiens")

        self.assertEqual(row_num, 0)  # only one data row in tsv file
        os.remove(tsv_filename)

        # Test json file of "homo_sapiens"
        json_filename = self.smash_path + "homo_sapiens/metadata_homo_sapiens.json"
        self.assertTrue(os.path.isfile(json_filename))
        with open(json_filename) as json_fp:
            species_metadada = json.load(json_fp)
        self.assertEqual(species_metadada["species"], "homo_sapiens")
        self.assertEqual(len(species_metadada["samples"]), 1)
        self.assertEqual(species_metadada["samples"][0]["refinebio_accession_code"], "GSM1361050")
        # os.remove(json_filename)

        # Test tsv file of "fake_species"
        tsv_filename = self.smash_path + "fake_species/metadata_fake_species.tsv"
        self.assertTrue(os.path.isfile(tsv_filename))
        with open(tsv_filename) as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            for row_num, row in enumerate(reader):
                self.assertEqual(row["refinebio_accession_code"], "E-GEOD-44719-GSM1089311")
                self.assertEqual(row["experiment_accession"], "E-GEOD-44719")
                self.assertEqual(
                    row["characteristic_cell population"], "IFNa DC"
                )  # ArrayExpress specific
                self.assertEqual(row["characteristic_dose"], "1 mL")  # ArrayExpress specific
                self.assertEqual(row["detection_percentage"], "98.44078")
                self.assertEqual(row["MetaSRA_age"], "3.0")

        self.assertEqual(row_num, 0)  # only one data row in tsv file
        os.remove(tsv_filename)

        # Test json file of "fake_species"
        json_filename = self.smash_path + "fake_species/metadata_fake_species.json"
        self.assertTrue(os.path.isfile(json_filename))

        with open(json_filename) as json_fp:
            species_metadada = json.load(json_fp)
        self.assertEqual(species_metadada["species"], "fake_species")
        self.assertEqual(len(species_metadada["samples"]), 1)
        self.assertEqual(
            species_metadada["samples"][0]["refinebio_accession_code"],
            "E-GEOD-44719-GSM1089311",
        )
        self.assertEqual(species_metadada["samples"][0]["MetaSRA_age"]["value"], 3.0)
        os.remove(json_filename)

    @tag("smasher")
    def test_bad_smash(self):
        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51081"
        experiment.save()

        result = ComputationalResult()
        result.save()

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        homo_sapiens.save()

        sample = Sample()
        sample.accession_code = "GSM1237810"
        sample.title = "GSM1237810"
        sample.organism = homo_sapiens
        sample.save()

        sample_annotation = SampleAnnotation()
        sample_annotation.data = {"hi": "friend"}
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
        computed_file.absolute_file_path = (
            "/home/user/data_store/BADSMASH/" + computed_file.filename
        )
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        sample = Sample()
        sample.accession_code = "GSM1237812"
        sample.title = "GSM1237812"
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
        computed_file.absolute_file_path = (
            "/home/user/data_store/BADSMASH/" + computed_file.filename
        )
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {"GSE51081": ["GSM1237810", "GSM1237812"]}
        ds.aggregate_by = "ALL"  # [ALL or SPECIES or EXPERIMENT]
        ds.scale_by = "NONE"  # [NONE or MINMAX or STANDARD or ROBUST]
        ds.email_address = "null@derp.com"
        # ds.email_address = "miserlou+heyo@gmail.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)
        self.assertTrue(ds.is_processed)

    @tag("smasher")
    def test_bad_overlap(self):
        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606)
        homo_sapiens.save()

        pj = ProcessorJob()
        pj.pipeline_applied = "SMASHER"
        pj.save()

        experiment = Experiment()
        experiment.accession_code = "GSE51081"
        experiment.save()

        result = ComputationalResult()
        result.save()

        # Now, make sure the bad can't zero this out.
        sample = Sample()
        sample.accession_code = "GSM999"
        sample.title = "GSM999"
        sample.organism = homo_sapiens
        sample.save()

        esa = ExperimentSampleAssociation()
        esa.experiment = experiment
        esa.sample = sample
        esa.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        computed_file = ComputedFile()
        computed_file.filename = "bad.PCL"
        computed_file.absolute_file_path = (
            "/home/user/data_store/BADSMASH/" + computed_file.filename
        )
        computed_file.result = result
        computed_file.size_in_bytes = 123
        computed_file.is_smashable = True
        computed_file.save()

        assoc = SampleComputedFileAssociation()
        assoc.sample = sample
        assoc.computed_file = computed_file
        assoc.save()

        ds = Dataset()
        ds.data = {"GSE51081": ["GSM1237810", "GSM1237812", "GSM999"]}
        ds.aggregate_by = "ALL"
        ds.scale_by = "NONE"
        ds.email_address = "null@derp.com"
        ds.quantile_normalize = False
        ds.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = pj
        pjda.dataset = ds
        pjda.save()

        final_context = smasher.smash(pj.pk, upload=False)
        ds = Dataset.objects.get(id=ds.id)

        self.assertEqual(len(final_context["final_frame"]), 4)

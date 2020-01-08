import hashlib
import os
import random
import shutil
import subprocess
from typing import Dict, List

from django.test import TestCase, tag

import numpy
import scipy.stats

from data_refinery_common.job_lookup import ProcessorEnum
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    OriginalFile,
    OriginalFileSampleAssociation,
    Pipeline,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    SurveyJob,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import salmon, tximport, utils


def prepare_organism_indices():
    c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

    # This is a lie, but this image doesn't have the dependencies for TRANSCRIPTOME_INDEX
    computational_result_short = ComputationalResult(processor=utils.find_processor("SALMON_QUANT"))
    computational_result_short.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = c_elegans
    organism_index.result = computational_result_short
    organism_index.absolute_directory_path = (
        "/home/user/data_store/salmon_tests/TRANSCRIPTOME_INDEX/SHORT"
    )
    organism_index.salmon_version = "salmon 0.13.1"
    organism_index.save()

    comp_file = ComputedFile()
    # This path will not be used because we already have the files extracted.
    comp_file.absolute_file_path = (
        "/home/user/data_store/salmon_tests/TRANSCRIPTOME_INDEX/SHORT/celgans_short.tar.gz"
    )
    comp_file.result = computational_result_short
    comp_file.size_in_bytes = 1337
    comp_file.sha1 = "ABC"
    comp_file.save()

    # This is a lie, but this image doesn't have the dependencies for TX_IMPORT
    computational_result_long = ComputationalResult(processor=utils.find_processor("SALMON_QUANT"))
    computational_result_long.save()


def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.save()

    c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

    samp = Sample()
    samp.accession_code = "SALMON"  # So the test files go to the right place
    samp.organism = c_elegans
    samp.source_database = "SRA"
    samp.technology = "RNA-SEQ"
    samp.save()

    prepare_organism_indices()

    og_file = OriginalFile()
    og_file.source_filename = "ERR1562482_1.fastq.gz"
    og_file.filename = "ERR1562482_1.fastq.gz"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR1562482_1.fastq.gz"
    og_file.is_downloaded = True
    og_file.save()

    og_file2 = OriginalFile()
    og_file2.source_filename = "ERR1562482_2.fastq.gz"
    og_file2.filename = "ERR1562482_2.fastq.gz"
    og_file2.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR1562482_2.fastq.gz"
    og_file2.is_downloaded = True
    og_file2.save()

    og_file_samp_assoc = OriginalFileSampleAssociation()
    og_file_samp_assoc.original_file = og_file
    og_file_samp_assoc.sample = samp
    og_file_samp_assoc.save()

    og_file_samp_assoc2 = OriginalFileSampleAssociation()
    og_file_samp_assoc2.original_file = og_file2
    og_file_samp_assoc2.sample = samp
    og_file_samp_assoc2.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file2
    assoc1.processor_job = pj
    assoc1.save()

    return pj, [og_file, og_file2]


def prepare_dotsra_job(filename="ERR1562482.sra"):
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.id = random.randint(111, 999999)
    pj.save()

    c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

    samp = Sample()
    samp.accession_code = "SALMON"  # So the test files go to the right place
    samp.organism = c_elegans
    samp.source_database = "SRA"
    samp.technology = "RNA-SEQ"
    samp.save()

    prepare_organism_indices()

    og_file = OriginalFile()
    og_file.source_filename = filename
    og_file.filename = filename
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/" + filename
    og_file.is_downloaded = True
    og_file.save()

    og_file_samp_assoc = OriginalFileSampleAssociation()
    og_file_samp_assoc.original_file = og_file
    og_file_samp_assoc.sample = samp
    og_file_samp_assoc.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj, [og_file]


def identical_checksum(filename1, filename2):
    """Confirm that the two files have identical checksum."""
    checksum_1 = hashlib.md5(open(filename1, "rb").read()).hexdigest()
    checksum_2 = hashlib.md5(open(filename2, "rb").read()).hexdigest()
    return checksum_1 == checksum_2


def strong_quant_correlation(ref_filename, output_filename):
    """Return true if both columns #3 and #4 (zero-indexed) of the two
    input quant files are strongly correlated (correlation >= 0.99).
    """
    ref_col34 = numpy.loadtxt(ref_filename, delimiter="\t", skiprows=1, usecols=(3, 4))
    ref_TPM = ref_col34[:, 0]
    ref_NumReads = ref_col34[:, 1]

    out_col34 = numpy.loadtxt(output_filename, delimiter="\t", skiprows=1, usecols=(3, 4))
    out_TPM = out_col34[:, 0]
    out_NumReads = out_col34[:, 1]

    TPM_stats = scipy.stats.spearmanr(ref_TPM, out_TPM)
    NumReads_stats = scipy.stats.spearmanr(ref_NumReads, out_NumReads)
    return TPM_stats.correlation >= 0.99 and NumReads_stats.correlation >= 0.99


class SalmonTestCase(TestCase):
    def setUp(self):
        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(
            name="CAENORHABDITIS_ELEGANS", taxonomy_id=6239, is_scientific_name=True
        )
        organism.save()

    @tag("salmon")
    def test_salmon(self):
        """Test the whole pipeline."""
        # Ensure any computed files from previous tests are removed.
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_job()
        job_context = salmon.salmon(job.pk)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertTrue(job.success)

        sample = files[0].samples.first()
        self.assertFalse(sample.is_processed)
        organism_index = job_context["quant_result"].organism_index
        self.assertEqual(organism_index.index_type, "TRANSCRIPTOME_SHORT")

    @tag("salmon")
    def test_no_salmon_on_geo(self):
        """Test that salmon won't be run on data coming from GEO."""
        # Ensure any computed files from previous tests are removed.
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_job()

        # We're expecting this processor job to fail, and when it does
        # it should clean up the original files that were for the
        # job. However we want to use these files in other tests, so
        # copy them so we can delete them without deleting the
        # originals.
        for original_file in OriginalFile.objects.all():
            new_path = original_file.absolute_file_path + "_copy"
            shutil.copyfile(original_file.absolute_file_path, new_path)
            original_file.absolute_file_path = new_path
            original_file.save()

        sample_object = Sample.objects.first()
        sample_object.source_database = "GEO"
        sample_object.save()

        job_context = salmon.salmon(job.pk)
        self.assertIsNotNone(job_context)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertFalse(job.success)
        self.assertEqual(
            job.failure_reason,
            ("The sample for this job either was not RNA-Seq or was not from the " "SRA database."),
        )
        self.assertTrue(job.no_retry)

        # Make sure the data got cleaned up, since the Janitor isn't
        # going to do it.
        for original_file in OriginalFile.objects.all():
            self.assertFalse(os.path.exists(original_file.absolute_file_path))

    @tag("salmon")
    def test_salmon_dotsra(self):
        """Test the whole pipeline."""
        # Ensure any computed files from previous tests are removed.
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_dotsra_job()
        job_context = salmon.salmon(job.pk)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertTrue(job.success)
        shutil.rmtree(job_context["work_dir"], ignore_errors=True)

    @tag("salmon")
    def test_salmon_dotsra_bad(self):
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_dotsra_job("i-dont-exist.sra")
        job_context = salmon.salmon(job.pk)
        self.assertIsNotNone(job_context)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertFalse(job.success)

    def check_salmon_quant(self, job_context, sample_dir):
        """Helper function that calls salmon._run_salmon and confirms
        strong correlation.
        """
        # Clean up if there were previous tests, but we still need that directory.
        shutil.rmtree(job_context["output_directory"], ignore_errors=True)
        os.makedirs(job_context["output_directory"], exist_ok=True)
        job_context = salmon._determine_index_length(job_context)
        job_context = salmon._find_or_download_index(job_context)

        job_context = salmon._run_salmon(job_context)
        job_context = salmon.get_tximport_inputs(job_context)
        job_context = salmon.tximport(job_context)
        output_quant_filename = os.path.join(job_context["output_directory"], "quant.sf")
        self.assertTrue(os.path.exists(output_quant_filename))

        # Confirm strong correlation between the new "quant.sf" and reference file
        ref_quant_filename = os.path.join(sample_dir, "ref_files/quant.sf")
        self.assertTrue(strong_quant_correlation(ref_quant_filename, output_quant_filename))

    @tag("salmon")
    def test_salmon_quant_one_sample_double_reads(self):
        """Test `salmon quant` on a sample that has double reads."""
        # Set up organism index database objects.
        prepare_organism_indices()

        # Create an Experiment that includes two samples.
        # (The first sample has test data available, but the second does not.)
        experiment_accession = "test_experiment"
        experiment = Experiment.objects.create(accession_code=experiment_accession)

        c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

        # test_sample record
        sample_accession = "test_sample"
        test_sample = Sample.objects.create(
            accession_code=sample_accession,
            organism=c_elegans,
            source_database="SRA",
            technology="RNA-SEQ",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=test_sample)
        # fake_sample record (created to prevent tximport step in this experiment)
        fake_sample = Sample.objects.create(
            accession_code="fake_sample", source_database="SRA", technology="RNA-SEQ"
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=fake_sample)

        experiment_dir = "/home/user/data_store/salmon_tests/test_experiment"

        og_read_1 = OriginalFile()
        og_read_1.absolute_file_path = os.path.join(experiment_dir, "raw/reads_1.fastq")
        og_read_1.filename = "reads_1.fastq"
        og_read_1.save()

        OriginalFileSampleAssociation.objects.create(
            original_file=og_read_1, sample=test_sample
        ).save()

        og_read_2 = OriginalFile()
        og_read_2.absolute_file_path = os.path.join(experiment_dir, "raw/reads_2.fastq")
        og_read_2.filename = "reads_1.fastq"
        og_read_2.save()

        OriginalFileSampleAssociation.objects.create(
            original_file=og_read_2, sample=test_sample
        ).save()

        sample_dir = os.path.join(experiment_dir, "test_sample")

        job_context = salmon._prepare_files(
            {
                "job_dir_prefix": "TEST",
                "job_id": "TEST",
                "job": ProcessorJob(),
                "pipeline": Pipeline(name="Salmon"),
                "computed_files": [],
                "original_files": [og_read_1, og_read_2],
            }
        )

        # Run salmon.
        self.check_salmon_quant(job_context, sample_dir)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample'.
        experiments_ready = salmon.get_tximport_inputs(job_context)["tximport_inputs"]
        self.assertEqual(len(experiments_ready), 0)

    @tag("salmon")
    def test_salmon_quant_two_samples_single_read(self):
        """Test `salmon quant` outputs on two samples that have single
        read and that belong to same experiment.
        """
        prepare_organism_indices()

        # Create one experiment and two related samples, based on:
        #   https://www.ncbi.nlm.nih.gov/sra/?term=SRP040623
        # (For testing purpose, only two of the four samples' data are included.)
        experiment_accession = "PRJNA242809"
        experiment = Experiment.objects.create(accession_code=experiment_accession)

        c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

        ## Sample 1
        sample1_accession = "SRR1206053"
        sample1 = Sample.objects.create(
            accession_code=sample1_accession,
            organism=c_elegans,
            source_database="SRA",
            technology="RNA-SEQ",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample1)

        experiment_dir = "/home/user/data_store/salmon_tests/PRJNA242809"

        og_file_1 = OriginalFile()
        og_file_1.absolute_file_path = os.path.join(experiment_dir, "raw/SRR1206053.fastq.gz")
        og_file_1.filename = "SRR1206053.fastq.gz"
        og_file_1.save()

        OriginalFileSampleAssociation.objects.create(original_file=og_file_1, sample=sample1).save()

        ## Sample 2
        sample2_accession = "SRR1206054"
        sample2 = Sample.objects.create(
            accession_code=sample2_accession,
            organism=c_elegans,
            source_database="SRA",
            technology="RNA-SEQ",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample2)

        og_file_2 = OriginalFile()
        og_file_2.absolute_file_path = os.path.join(experiment_dir, "raw/SRR1206054.fastq.gz")
        og_file_2.filename = "SRR1206054.fastq.gz"
        og_file_2.save()

        OriginalFileSampleAssociation.objects.create(original_file=og_file_2, sample=sample2).save()

        # Test `salmon quant` on sample1 (SRR1206053)
        sample1_dir = os.path.join(experiment_dir, sample1_accession)

        job1_context = salmon._prepare_files(
            {
                "job_dir_prefix": "TEST",
                "job_id": "TEST",
                "pipeline": Pipeline(name="Salmon"),
                "computed_files": [],
                "original_files": [og_file_1],
            }
        )

        # Check quant.sf in `salmon quant` output dir of sample1
        self.check_salmon_quant(job1_context, sample1_dir)
        # Confirm that this experiment is not ready for tximport yet.
        experiments_ready = salmon.get_tximport_inputs(job1_context)["tximport_inputs"]
        self.assertEqual(len(experiments_ready), 0)
        # This job should not have produced any tximport output
        # because the other sample isn't ready yet.
        self.assertFalse(os.path.exists(os.path.join(job1_context["work_dir"], "txi_out.RDS")))

        # Now run `salmon quant` on sample2 (SRR1206054) too
        sample2_dir = os.path.join(experiment_dir, sample2_accession)
        job2_context = salmon._prepare_files(
            {
                "job_dir_prefix": "TEST2",
                "job_id": "TEST2",
                "pipeline": Pipeline(name="Salmon"),
                "computed_files": [],
                "original_files": [og_file_2],
            }
        )

        # Clean up tximport output:
        rds_filename = os.path.join(job2_context["work_dir"], "txi_out.RDS")
        if os.path.isfile(rds_filename):
            os.remove(rds_filename)

        # Check quant.sf in `salmon quant` output dir of sample2
        self.check_salmon_quant(job2_context, sample2_dir)

        # rds_filename should have been generated by tximport at this point.
        # Note: `tximport` step is launched by subprocess module in Python.
        # If input "quant.sf" files are too large, we may have to wait for
        # a few seconds before testing the existence of rds_filename.
        self.assertTrue(os.path.exists(rds_filename))

        for computed_file in job2_context["computed_files"]:
            if computed_file.filename[-4:] == ".RDS":
                rds_file_path = computed_file.absolute_file_path

        cmd_tokens = [
            "/usr/bin/Rscript",
            "--vanilla",
            "/home/user/data_refinery_workers/processors/test_tximport.R",
            "--txi_out",
            rds_file_path,
            "--gene2txmap",
            job2_context["genes_to_transcripts_path"],
        ]

        tximport_test_result = subprocess.run(
            cmd_tokens, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        if tximport_test_result.returncode != 0:
            # If the exit code is not 0 then tximport failed so fail the tests.
            self.assertTrue(False)

        # Check the individual files
        self.assertTrue(len(job2_context["individual_files"]), 2)
        for file in job2_context["individual_files"]:
            self.assertTrue(os.path.isfile(file.absolute_file_path))

    @tag("salmon")
    def test_get_tximport_inputs(self):
        """"Tests that tximport only considers RNA-Seq samples from GEO.
        """
        # Create one experiment and two related samples, based on:
        #   https://www.ncbi.nlm.nih.gov/sra/?term=SRP040623
        # (We don't need any original files because
        # get_tximport_inputs doesn't consider them.)
        experiment_accession = "PRJNA242809"
        experiment = Experiment.objects.create(accession_code=experiment_accession)

        c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

        ## Sample 1
        sample1_accession = "SRR1206053"
        sample1 = Sample.objects.create(accession_code=sample1_accession, organism=c_elegans)
        sample1.source_database = "GEO"
        sample1.technology = "RNA-SEQ"
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample1)

        ## Sample 2
        sample2_accession = "SRR1206054"
        sample2 = Sample.objects.create(accession_code=sample2_accession, organism=c_elegans)
        sample2.source_database = "GEO"
        sample2.technology = "RNA-SEQ"
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample2)

        computational_result1 = ComputationalResult(processor=utils.find_processor("SALMON_QUANT"))
        computational_result1.save()

        sample_result_assoc = SampleResultAssociation(sample=sample1, result=computational_result1)
        sample_result_assoc.save()

        comp_file = ComputedFile()
        comp_file.absolute_file_path = "/doesnt/matter"
        comp_file.result = computational_result1
        comp_file.size_in_bytes = 1337
        comp_file.sha1 = "ABC"
        comp_file.s3_key = "key"
        comp_file.s3_bucket = "bucket"
        comp_file.save()

        computational_result2 = ComputationalResult(processor=utils.find_processor("SALMON_QUANT"))
        computational_result2.save()

        sample_result_assoc = SampleResultAssociation(sample=sample2, result=computational_result2)
        sample_result_assoc.save()

        comp_file = ComputedFile()
        comp_file.absolute_file_path = "/doesnt/matter"
        comp_file.result = computational_result2
        comp_file.size_in_bytes = 1337
        comp_file.sha1 = "ABC"
        comp_file.s3_key = "key"
        comp_file.s3_bucket = "bucket"
        comp_file.save()

        quantified_experiments = salmon.get_tximport_inputs({"sample": sample1})["tximport_inputs"]

        self.assertEqual({}, quantified_experiments)


class SalmonToolsTestCase(TestCase):
    """Test SalmonTools command."""

    def setUp(self):
        self.test_dir = "/home/user/data_store/salmontools/"
        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(
            name="CAENORHABDITIS_ELEGANS", taxonomy_id=6239, is_scientific_name=True
        )
        organism.save()

    @tag("salmon")
    def test_double_reads(self):
        """Test outputs when the sample has both left and right reads."""
        job_context = {
            "job_id": 123,
            "job": ProcessorJob(),
            "pipeline": Pipeline(name="Salmon"),
            "input_file_path": self.test_dir + "double_input/reads_1.fastq",
            "input_file_path_2": self.test_dir + "double_input/reads_2.fastq",
            "salmontools_directory": self.test_dir + "double_salmontools/",
            "salmontools_archive": self.test_dir + "salmontools-result.tar.gz",
            "output_directory": self.test_dir + "double_output/",
            "computed_files": [],
        }
        os.makedirs(job_context["salmontools_directory"], exist_ok=True)

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        sample = Sample()
        sample.organism = homo_sapiens
        sample.save()
        job_context["sample"] = sample

        salmon._run_salmontools(job_context)

        # Confirm job status
        self.assertTrue(job_context["success"])

        # Unpack result for checking
        os.system("gunzip " + job_context["salmontools_directory"] + "*.gz")

        # Check two output files
        output_file1 = job_context["salmontools_directory"] + "unmapped_by_salmon_1.fa"
        expected_output_file1 = self.test_dir + "expected_double_output/unmapped_by_salmon_1.fa"
        self.assertTrue(identical_checksum(output_file1, expected_output_file1))

        output_file2 = job_context["salmontools_directory"] + "unmapped_by_salmon_2.fa"
        expected_output_file2 = self.test_dir + "expected_double_output/unmapped_by_salmon_2.fa"
        self.assertTrue(identical_checksum(output_file2, expected_output_file2))

    @tag("salmon")
    def test_single_read(self):
        """Test outputs when the sample has one read only."""
        job_context = {
            "job_id": 456,
            "job": ProcessorJob(),
            "pipeline": Pipeline(name="Salmon"),
            "input_file_path": self.test_dir + "single_input/single_read.fastq",
            "output_directory": self.test_dir + "single_output/",
            "salmontools_directory": self.test_dir + "single_salmontools/",
            "salmontools_archive": self.test_dir + "salmontools-result.tar.gz",
            "computed_files": [],
        }
        os.makedirs(job_context["salmontools_directory"], exist_ok=True)

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        sample = Sample()
        sample.organism = homo_sapiens
        sample.save()
        job_context["sample"] = sample

        salmon._run_salmontools(job_context)

        # Confirm job status
        self.assertTrue(job_context["success"])

        # Unpack result for checking
        os.system("gunzip " + job_context["salmontools_directory"] + "*.gz")

        # Check output file
        output_file = job_context["salmontools_directory"] + "unmapped_by_salmon.fa"
        expected_output_file = self.test_dir + "expected_single_output/unmapped_by_salmon.fa"
        self.assertTrue(identical_checksum(output_file, expected_output_file))


class DetermineIndexLengthTestCase(TestCase):
    """Test salmon._determine_index_length function, which gets the salmon index length of a sample.
    For now, these tests only ensure that the output of the new faster salmon index function match
    that of the old one for the test data."""

    def setUp(self):
        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(
            name="CAENORHABDITIS_ELEGANS", taxonomy_id=6239, is_scientific_name=True
        )
        organism.save()

    @tag("salmon")
    def test_salmon_determine_index_length_single_read(self):
        """Test that the right length is calculated when the sample has one read."""
        job, files = prepare_job()

        job_context = salmon._set_job_prefix(
            {"original_files": [files[0]], "job_id": job.id, "job": job}
        )
        job_context = salmon._prepare_files(job_context)
        results = salmon._determine_index_length(job_context)

        self.assertEqual(results["index_length_raw"], 41)
        self.assertEqual(results["index_length"], "short")

    @tag("salmon")
    def test_salmon_determine_index_length_double_read(self):
        """Test that the right length is calculated when the sample has two reads."""
        job, files = prepare_job()

        job_context = salmon._set_job_prefix(
            {"original_files": files, "job_id": job.id, "job": job}
        )
        job_context = salmon._prepare_files(job_context)
        results = salmon._determine_index_length(job_context)

        self.assertEqual(results["index_length_raw"], 41)
        self.assertEqual(results["index_length"], "short")


class RuntimeProcessorTest(TestCase):
    """Test the four processors hosted inside "Salmon" docker container."""

    @tag("salmon")
    def test_tximport(self):
        self.assertEqual(Processor.objects.count(), 0)  # No processor yet

        proc_key = "TXIMPORT"
        tximport_processor = utils.find_processor(proc_key)
        self.assertEqual(Processor.objects.count(), 1)  # New processor created

        # Validate some information of the new processor
        self.assertEqual(tximport_processor.name, ProcessorEnum[proc_key].value["name"])
        self.assertEqual(tximport_processor.version, get_env_variable("SYSTEM_VERSION"))
        self.assertEqual(
            tximport_processor.docker_image, ProcessorEnum[proc_key].value["docker_img"]
        )
        self.assertEqual(tximport_processor.environment["os_distribution"], utils.get_os_distro())

        os_pkg_name = "r-base"
        self.assertEqual(
            tximport_processor.environment["os_pkg"][os_pkg_name],
            utils.get_os_pkgs([os_pkg_name])[os_pkg_name],
        )

        pip_pkg_name = "data-refinery-common"
        self.assertEqual(
            tximport_processor.environment["python"][pip_pkg_name],
            utils.get_pip_pkgs([pip_pkg_name])[pip_pkg_name],
        )

        r_pkg_names = ["Bioconductor", "tximport"]
        r_pkg_info = utils.get_r_pkgs(r_pkg_names)
        for r_pkg in r_pkg_names:
            self.assertEqual(tximport_processor.environment["R"][r_pkg], r_pkg_info[r_pkg])

        # Confirm that there is only one processor in one runtime environment
        for i in range(3):
            proc2 = utils.find_processor(proc_key)
            self.assertEqual(Processor.objects.count(), 1)  # No new processor
            self.assertEqual(tximport_processor, proc2)  # Same processor instance

    @tag("salmon")
    def test_salmon_quant(self):
        self.assertEqual(Processor.objects.count(), 0)  # No processor yet
        proc_key = "SALMON_QUANT"
        sq_processor = utils.find_processor(proc_key)
        self.assertEqual(Processor.objects.count(), 1)  # New processor created

        # Validate some information of the new processor
        self.assertEqual(sq_processor.name, ProcessorEnum[proc_key].value["name"])

        cmd_str = "salmon --version"
        cmd_output = utils.get_cmd_lines([cmd_str])[cmd_str]
        self.assertEqual(sq_processor.environment["cmd_line"][cmd_str], cmd_output)

    @tag("salmon")
    def test_salmontools(self):
        self.assertEqual(Processor.objects.count(), 0)  # No processor yet
        proc_key = "SALMONTOOLS"
        st_processor = utils.find_processor(proc_key)
        self.assertEqual(Processor.objects.count(), 1)  # New processor created

        # Validate some information of the new processor
        self.assertEqual(st_processor.name, ProcessorEnum[proc_key].value["name"])

        cmd_str = "salmontools --version"
        cmd_output = utils.get_cmd_lines([cmd_str])[cmd_str]
        self.assertEqual(st_processor.environment["cmd_line"][cmd_str], cmd_output)

    @tag("salmon")
    def test_exception_handler(self):
        """Test utils.handle_processor_expcetion to confirm that the
        exception is handled correctly.
        """
        proc_key = "foobar"  # This processor key does NOT exist!
        job_context = dict()
        job_context["job"] = ProcessorJob.objects.create()
        job_context["success"] = True
        try:
            utils.find_processor(proc_key)
        except Exception as e:
            utils.handle_processor_exception(job_context, proc_key, e)

        # Failed job because "foobar" is an invalid processor key
        self.assertEqual(job_context["success"], False)
        self.assertEqual(job_context["job"].failure_reason, "Failed to set processor: 'foobar'")

    @tag("salmon")
    def test_salmontools_with_bad_processor(self):
        """Test salmontools with a bad processor key."""
        test_dir = "/home/user/data_store/salmontools/"
        job_context = {
            "job_id": 123,
            "job": ProcessorJob.objects.create(),
            "pipeline": Pipeline(name="Salmon"),
            "input_file_path": test_dir + "double_input/reads_1.fastq",
            "input_file_path_2": test_dir + "double_input/reads_2.fastq",
            "salmontools_directory": test_dir + "double_salmontools/",
            "salmontools_archive": test_dir + "salmontools-result.tar.gz",
            "output_directory": test_dir + "double_output/",
        }
        os.makedirs(job_context["salmontools_directory"], exist_ok=True)
        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        sample = Sample()
        sample.organism = homo_sapiens
        sample.save()
        job_context["sample"] = sample

        # Set the wrong yml filename on purpose to mess up Salmontools processor
        original_yml_file = ProcessorEnum["SALMONTOOLS"].value["yml_file"]
        ProcessorEnum["SALMONTOOLS"].value["yml_file"] = "foobar.yml"

        salmon._run_salmontools(job_context)
        self.assertEqual(job_context["success"], False)
        self.assertTrue(job_context["job"].failure_reason.startswith("Failed to set processor:"))

        # Change yml filename back
        ProcessorEnum["SALMONTOOLS"].value["yml_file"] = original_yml_file


def create_tximport_job_context(
    complete_accessions: List[str], incomplete_accessions: List[str], salmon_version="salmon 0.13.1"
) -> Dict:
    """Create an experiment and associated objects and run tximport on it.

    Creates a sample for each accession contained in either input
    list. The samples in complete_accessions will be simlulated as
    already having salmon quant run on them. The samples in
    incomplete_accessions won't.
    """
    # Create the experiment
    experiment_accession = "SRP095529"
    data_dir = "/home/user/data_store/"
    experiment_dir = data_dir + experiment_accession
    experiment = Experiment.objects.create(accession_code=experiment_accession)

    zebrafish = Organism.get_object_for_name("DANIO_RERIO")

    # This is a lie, but this image doesn't have the dependencies for TRANSCRIPTOME_INDEX
    computational_result_short = ComputationalResult(processor=utils.find_processor("SALMON_QUANT"))
    computational_result_short.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = zebrafish
    organism_index.result = computational_result_short
    organism_index.absolute_directory_path = "/home/user/data_store/ZEBRAFISH_INDEX/SHORT"
    organism_index.salmon_version = salmon_version
    organism_index.save()

    comp_file = ComputedFile()
    # This path will not be used because we already have the files extracted.
    comp_file.absolute_file_path = (
        "/home/user/data_store/ZEBRAFISH_INDEX/SHORT/zebrafish_short.tar.gz"
    )
    comp_file.result = computational_result_short
    comp_file.size_in_bytes = 1337
    comp_file.sha1 = "ABC"
    comp_file.s3_key = "key"
    comp_file.s3_bucket = "bucket"
    comp_file.save()

    for accession_code in incomplete_accessions:
        last_sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database="SRA",
            technology="RNA-SEQ",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=last_sample)

    # Create tximport result and files
    quant_processor = utils.find_processor("SALMON_QUANT")
    utils.find_processor("TXIMPORT")

    # Create the already processed samples along with their
    # ComputationalResults and ComputedFiles. They don't need
    # original files for this test because we aren't going to run
    # salmon quant on them.
    for accession_code in complete_accessions:
        sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database="SRA",
            technology="RNA-SEQ",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

        if accession_code == "SRR5125622":
            current_sample = sample

        # Create and associate quant result and files.
        quant_result = ComputationalResult()
        quant_result.is_ccdl = True
        quant_result.processor = quant_processor
        quant_result.organism_index = organism_index
        quant_result.save()

        kv = ComputationalResultAnnotation()
        kv.data = {"index_length": "short"}
        kv.result = quant_result
        kv.is_public = True
        kv.save()

        # In prod the filename pattern will involve the timestamp
        # but here we're using the accession code so we can find
        # the archive file for the current sample.
        archive_filename = "result-" + accession_code + ".tar.gz"
        archive_file = ComputedFile()
        archive_file.filename = archive_filename
        archive_file.absolute_file_path = os.path.join(experiment_dir, archive_filename)
        archive_file.is_public = False
        archive_file.is_smashable = False
        archive_file.is_qc = False
        archive_file.result = quant_result
        archive_file.size_in_bytes = 12345
        archive_file.save()

        quant_file = ComputedFile()
        quant_file.filename = "quant.sf"
        quant_file.absolute_file_path = (
            experiment_dir + "/quant_files/" + accession_code + "_output/quant.sf"
        )
        quant_file.is_public = False
        quant_file.is_smashable = False
        quant_file.is_qc = False
        quant_file.result = quant_result
        quant_file.size_in_bytes = 12345
        quant_file.s3_bucket = "bucket"
        quant_file.s3_key = "key"
        quant_file.save()

        SampleResultAssociation.objects.get_or_create(sample=sample, result=quant_result)

    # Processor jobs need at least one original file associated with
    # them so they know what they're processing.
    current_og = OriginalFile()
    current_og.absolute_file_path = os.path.join(experiment_dir, "SRR5125622.fastq.gz")
    current_og.filename = "SRR5125622.fastq.gz"
    current_og.save()

    OriginalFileSampleAssociation.objects.create(
        original_file=current_og, sample=current_sample
    ).save()

    pj = ProcessorJob()
    pj.pipeline_applied = "TXIMPORT"
    pj.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = current_og
    assoc1.processor_job = pj
    assoc1.save()

    # Prep our job context
    job_context = tximport._prepare_files(
        {
            "job_dir_prefix": "TEST3",
            "job_id": "TEST3",
            "job": pj,
            "index_directory": organism_index.absolute_directory_path,
            "pipeline": Pipeline(name="Salmon"),
            "computed_files": [],
            "original_files": [current_og],
        }
    )

    # We don't have the raw file to run _determine_index_length so
    # just pick one, it doesn't matter that much because we aren't
    # checking the output data.
    job_context["index_length"] = "short"

    return job_context


def run_tximport_for_job_context(job_context: Dict) -> Dict:
    job_context = salmon._find_or_download_index(job_context)

    job_context = salmon.get_tximport_inputs(job_context)
    job_context = salmon.tximport(job_context)
    job_context = utils.end_job(job_context)

    return job_context


def run_tximport_at_progress_point(
    complete_accessions: List[str], incomplete_accessions: List[str], salmon_version="salmon 0.13.1"
) -> Dict:
    job_context = create_tximport_job_context(
        complete_accessions, incomplete_accessions, salmon_version
    )

    return run_tximport_for_job_context(job_context)


class EarlyTximportTestCase(TestCase):
    """Tests that running tximport early works when it should be triggered.

    Some experiments are going to have samples that can't be
    processed. This means that tximport can't run on those unless we
    tell it to. We have a tximport only job to do this, but it should
    only be run on experiments with at least 25 samples where at least
    80% of the samples have been processed. We therefore run the
    tximport job on an experiment that is ready for it, one that has
    too few samples, and one that has too low of a copmpletion
    percent.
    """

    @tag("salmon")
    def test_early_tximport(self):
        """Tests that running tximport early works correctly.

        Makes sure that when we should in fact run tximport early that
        we do so, it works, and that it works even if there already is
        an existing tximport result for the experiment.

        So we run tximport on an experiment with 5 samples that aren't
        yet processed and 20 samples that are. By having 20/25 samples
        complete, we're just past both the numerical and percent
        cutoffs.
        """
        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        self.assertTrue("tximported" in job_context)

        rds_file = ComputedFile.objects.get(filename="txi_out.RDS")

        for accession_code in complete_accessions:
            # Check to make sure that all the associations were madep
            # correctly. These queries will fail if they weren't.
            tpm_file = ComputedFile.objects.get(
                filename=accession_code + "_output_gene_lengthScaledTPM.tsv"
            )
            sample = tpm_file.samples.first()
            SampleComputedFileAssociation.objects.get(sample=sample, computed_file=tpm_file)
            SampleComputedFileAssociation.objects.get(sample=sample, computed_file=rds_file)
            self.assertTrue(sample.is_processed)

        # Make sure that these samples actually were ignored.
        for accession_code in incomplete_accessions:
            sample = Sample.objects.get(accession_code=accession_code)
            self.assertEqual(sample.computed_files.count(), 0)

    @tag("salmon")
    def test_version_filter(self):
        """Tests that we don't run tximport on old salmon versions.
        """
        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = run_tximport_at_progress_point(
            complete_accessions, incomplete_accessions, salmon_version="salmon 0.9.1"
        )

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertFalse("tximported" in job_context)

    @tag("salmon")
    def test_tximport_percent_cutoff(self):
        """Tests logic for determining if tximport should be run early.

        This test is verifying that tximport won't run if the
        experiment has too few samples in it to be elegible.

        So we run tximport on an experiment with 6 samples that aren't
        yet processed and 20 samples that are. By having 20/26 samples
        complete, we're past the numerical cutoff, but not past the
        percent cutoff.
        """
        # Accessions SRR5125615-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125615",
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertFalse("tximported" in job_context)

        # Make sure that these samples actually were ignored.
        for accession_code in incomplete_accessions:
            sample = Sample.objects.get(accession_code=accession_code)
            self.assertEqual(sample.computed_files.count(), 0)

    @tag("salmon")
    def test_tximport_numerical_cutoff(self):
        """Tests logic for determining if tximport should be run early.

        This test is verifying that tximport won't run if the
        experiment has too few samples in it to be elegible.

        So we run tximport on an experiment with 5 samples that aren't
        yet processed and 19 that are. By having 19/24 samples
        complete, we're past the percent cutoff, but not past the
        numerical cutoff.
        """
        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertFalse("tximported" in job_context)

        # Make sure that these samples actually were ignored.
        for accession_code in incomplete_accessions:
            sample = Sample.objects.get(accession_code=accession_code)
            self.assertEqual(sample.computed_files.count(), 0)

    @tag("salmon")
    def test_missing_computed_file(self):
        """Tests that tximport will ignore computed files that were not uploaded to S3.

        We will do this by setting the s3_bucket and s3_key fields to
        None for one of the complete_accessions so that it will get
        filtered out and then we won't make the threshold for
        tximport.
        """
        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = create_tximport_job_context(complete_accessions, incomplete_accessions)

        result = Sample.objects.filter(accession_code="SRR5125621").first().results.first()
        computed_file = ComputedFile.objects.get(result=result, filename="quant.sf")
        computed_file.s3_bucket = None
        computed_file.s3_key = None
        computed_file.save()

        job_context = run_tximport_for_job_context(job_context)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertFalse("tximported" in job_context)

        # Make sure that these samples actually were ignored.
        for accession_code in incomplete_accessions:
            sample = Sample.objects.get(accession_code=accession_code)
            self.assertEqual(sample.computed_files.count(), 0)

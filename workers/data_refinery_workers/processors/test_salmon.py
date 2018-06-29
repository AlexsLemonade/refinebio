import hashlib
import os
import shutil
import numpy
import scipy.stats
from contextlib import closing
from django.test import TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    Organism,
    OrganismIndex,
    ComputedFile,
    ComputationalResult,
    Processor,
    Pipeline,
    Sample,
    Experiment,
    ExperimentSampleAssociation,
    ProcessorJobOriginalFileAssociation,
    OriginalFileSampleAssociation
)
from data_refinery_workers.processors import salmon, utils
from data_refinery_workers._version import __version__


def setUpModule():
    for program in ["Salmon Quant", "Tximport", "MultiQC", "Salmontools"]:
        processor_name = program + " " + __version__
        Processor.objects.create(name=processor_name)


def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.save()

    c_elegans = Organism.get_object_for_name("CAENORHABDITIS_ELEGANS")

    samp = Sample()
    samp.accession_code = "SALMON" # So the test files go to the right place
    samp.organism = c_elegans
    samp.save()

    computational_result = ComputationalResult(processor=Processor.objects.first())
    computational_result.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = c_elegans
    organism_index.result = computational_result
    organism_index.save()

    comp_file = ComputedFile()
    comp_file.absolute_file_path = "/home/user/data_store/processed/TEST/TRANSCRIPTOME_INDEX/Caenorhabditis_elegans_short_1527089586.tar.gz"
    comp_file.result = computational_result
    comp_file.calculate_size()
    comp_file.calculate_sha1()
    comp_file.save()

    og_file = OriginalFile()
    og_file.source_filename = "ERR1562482_1.fastq.gz"
    og_file.filename = "ERR1562482_1.fastq.gz"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR1562482_1.fastq.gz"
    og_file.save()

    og_file2 = OriginalFile()
    og_file2.source_filename = "ERR1562482_2.fastq.gz"
    og_file2.filename = "ERR1562482_2.fastq.gz"
    og_file2.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR1562482_2.fastq.gz"
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


def identical_checksum(filename1, filename2):
    """Confirm that the two files have identical checksum."""
    checksum_1 = hashlib.md5(open(filename1, 'rb').read()).hexdigest()
    checksum_2 = hashlib.md5(open(filename2, 'rb').read()).hexdigest()
    return checksum_1 == checksum_2


def strong_quant_correlation(ref_filename, output_filename):
    """Return true if both columns #3 and #4 (zero-indexed) of the two
    input quant files are strongly correlated (correlation >= 0.99).
    """
    ref_col34 = numpy.loadtxt(ref_filename, delimiter='\t', skiprows=1, usecols=(3,4))
    ref_TPM = ref_col34[:, 0]
    ref_NumReads = ref_col34[:, 1]

    out_col34 = numpy.loadtxt(output_filename, delimiter='\t', skiprows=1, usecols=(3,4))
    out_TPM = out_col34[:, 0]
    out_NumReads = out_col34[:, 1]

    TPM_stats = scipy.stats.spearmanr(ref_TPM, out_TPM)
    NumReads_stats = scipy.stats.spearmanr(ref_NumReads, out_NumReads)
    return TPM_stats.correlation >= 0.99 and NumReads_stats.correlation >= 0.99


class SalmonTestCase(TestCase):
    def setUp(self):
        self.test_dir = '/home/user/data_store/salmon_tests'

    @tag('salmon')
    def test_salmon(self):
        """Test the whole pipeline."""
        # Ensure any computed files from previous tests are removed.
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_job()
        salmon.salmon(job.pk)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertTrue(job.success)

    def chk_salmon_quant(self, job_context, sample_dir):
        """Helper function that calls salmon._run_salmon and confirms
        strong correlation.
        """

        shutil.rmtree(job_context['output_directory'], ignore_errors=True)  # clean up
        salmon._run_salmon(job_context, skip_processed=False)
        output_quant_filename = os.path.join(job_context['output_directory'], 'quant.sf')
        self.assertTrue(os.path.exists(output_quant_filename))

        # Confirm strong correlation between the new "quant.sf" and reference file
        ref_quant_filename = os.path.join(sample_dir, 'ref_files/quant.sf')
        self.assertTrue(strong_quant_correlation(ref_quant_filename, output_quant_filename))

    @tag('salmon')
    def test_salmon_quant_one_sample_double_reads(self):
        """Test `salmon quant` on a sample that has double reads."""

        # Create an Experiment that includes two samples.
        # (The first sample has test data available, but the second does not.)
        experiment_accession = 'test_experiment'
        experiment = Experiment.objects.create(accession_code=experiment_accession)
        # test_sample record
        sample_accession = 'test_sample'
        test_sample = Sample.objects.create(accession_code=sample_accession)
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=test_sample)
        # fake_sample record (created to prevent tximport step in this experiment)
        fake_sample = Sample.objects.create(accession_code='fake_sample')
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=fake_sample)

        experiment_dir = os.path.join(self.test_dir, experiment_accession)
        sample_dir = os.path.join(experiment_dir, 'test_sample')

        job_context = {
            'job_id': 1,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'sample': test_sample,
            'input_file_path': os.path.join(experiment_dir, 'raw/reads_1.fastq'),
            'input_file_path_2': os.path.join(experiment_dir, 'raw/reads_2.fastq'),
            'index_directory': os.path.join(experiment_dir, 'index'),
            'output_directory': os.path.join(sample_dir, 'processed'),
            'success': True
        }
        # Check quant.sf in `salmon quant` output dir
        self.chk_salmon_quant(job_context, sample_dir)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample'.
        experiments_ready = salmon._get_salmon_completed_exp_dirs(job_context)
        self.assertEqual(len(experiments_ready), 0)

    @tag('salmon')
    def test_salmon_quant_two_samples_single_read(self):
        """Test `salmon quant` outputs on two samples that have single
        read and that belong to same experiment.
        """
        # Create one experiment and two related samples, based on:
        #   https://www.ncbi.nlm.nih.gov/sra/?term=SRP040623
        # (For testing purpose, only two of the four samples' data are included.)
        experiment_accession = 'PRJNA242809'
        experiment = Experiment.objects.create(accession_code=experiment_accession)

        sample1_accession = 'SRR1206053'
        sample1 = Sample.objects.create(accession_code=sample1_accession)
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample1)

        sample2_accession = 'SRR1206054'
        sample2 = Sample.objects.create(accession_code=sample2_accession)
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample2)

        experiment_dir = os.path.join(self.test_dir, experiment_accession)
        # Clean up tximport output:
        rds_filename = os.path.join(experiment_dir, 'txi_out.RDS')
        if (os.path.isfile(rds_filename)):
            os.remove(rds_filename)

        # Test `salmon quant` on sample1 (SRR1206053)
        sample1_dir = os.path.join(experiment_dir, sample1_accession)
        job1_context = {
            'job_id': 1,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'sample': sample1,
            'input_file_path': os.path.join(experiment_dir, 'raw/SRR1206053.fastq.gz'),
            'index_directory': os.path.join(experiment_dir, 'index'),
            'genes_to_transcripts_path': os.path.join(experiment_dir, 'index',
                                                      'genes_to_transcripts.txt'),
            'output_directory': os.path.join(sample1_dir, 'processed'),
            'success': True
        }

        # Check quant.sf in `salmon quant` output dir of sample1
        self.chk_salmon_quant(job1_context, sample1_dir)
        # Confirm that this experiment is not ready for tximport yet.
        experiments_ready = salmon._get_salmon_completed_exp_dirs(job1_context)
        self.assertEqual(len(experiments_ready), 0)
        self.assertFalse(os.path.exists(rds_filename))

        # Now run `salmon quant` on sample2 (SRR1206054) too
        sample2_dir = os.path.join(experiment_dir, sample2_accession)
        job2_context = {
            'job_id': 2,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'sample': sample2,
            'input_file_path': os.path.join(experiment_dir, 'raw/SRR1206054.fastq.gz'),
            'index_directory': os.path.join(experiment_dir, 'index'),
            'genes_to_transcripts_path': os.path.join(experiment_dir, 'index',
                                                      'genes_to_transcripts.txt'),
            'output_directory': os.path.join(sample2_dir, 'processed'),
            'success': True
        }

        # Check quant.sf in `salmon quant` output dir of sample2
        self.chk_salmon_quant(job2_context, sample2_dir)
        # Confirm that this experiment is ready for tximport now:
        experiments_ready = salmon._get_salmon_completed_exp_dirs(job2_context)
        self.assertEqual(len(experiments_ready), 1)
        self.assertEqual(experiments_ready[0], experiment_dir)

        # rds_filename should have been generated bytximport at this point.
        # Note: `tximport` step is launched by subprocess module in Python.
        # If input "quant.sf" files are too large, we may have to wait for
        # a few seconds before testing the existence of rds_filename.
        self.assertTrue(os.path.exists(rds_filename))

    def test_fastqc(self):

        job, og_files = prepare_job()
        win_context = {
            'job': job,
            'job_id': 789,
            'pipeline': Pipeline(name="Salmon"),
            'qc_directory': "/home/user/data_store/raw/TEST/SALMON/qc",
            'original_files': og_files,
            'success': True
        }

        # Ensure clean testdir
        shutil.rmtree(win_context['qc_directory'], ignore_errors=True)
        os.makedirs(win_context['qc_directory'], exist_ok=True)
        win_context = salmon._prepare_files(win_context)

        win = salmon._run_fastqc(win_context)
        self.assertTrue(win['success'])
        win = salmon._run_multiqc(win_context)
        self.assertTrue(win['success'])

        for file in win['qc_files']:
            self.assertTrue(os.path.isfile(file.absolute_file_path))

        fail_context = {
            'job': job,
            'job_id': 'hippityhoppity',
            'pipeline': Pipeline(name="Salmon"),
            'qc_directory': "/home/user/data_store/raw/TEST/SALMON/derp",
            'original_files': [],
            'success': True
        }
        fail = salmon._run_fastqc(fail_context)
        self.assertFalse(fail['success'])


class SalmonToolsTestCase(TestCase):
    """Test SalmonTools command."""

    def setUp(self):
        self.test_dir = '/home/user/data_store/salmontools/'

    @tag('salmon')
    def test_double_reads(self):
        """Test outputs when the sample has both left and right reads."""
        job_context = {
            'job_id': 123,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'input_file_path': self.test_dir + 'double_input/reads_1.fastq',
            'input_file_path_2': self.test_dir + 'double_input/reads_2.fastq',
            'output_directory': self.test_dir + 'double_output/'
        }

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        sample = Sample()
        sample.organism = homo_sapiens
        sample.save()
        job_context["sample"] = sample

        salmon._run_salmontools(job_context, False)

        # Confirm job status
        self.assertTrue(job_context["success"])

        # Check two output files
        output_file1 = self.test_dir + 'double_output/unmapped_by_salmon_1.fa'
        expected_output_file1 = self.test_dir + 'expected_double_output/unmapped_by_salmon_1.fa'
        self.assertTrue(identical_checksum(output_file1, expected_output_file1))

        output_file2 = self.test_dir + 'double_output/unmapped_by_salmon_2.fa'
        expected_output_file2 = self.test_dir + 'expected_double_output/unmapped_by_salmon_2.fa'
        self.assertTrue(identical_checksum(output_file2, expected_output_file2))

    @tag('salmon')
    def test_single_read(self):
        """Test outputs when the sample has one read only."""
        job_context = {
            'job_id': 456,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'input_file_path': self.test_dir + 'single_input/single_read.fastq',
            'output_directory': self.test_dir + 'single_output/'
        }

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        sample = Sample()
        sample.organism = homo_sapiens
        sample.save()
        job_context["sample"] = sample

        salmon._run_salmontools(job_context, False)

        # Confirm job status
        self.assertTrue(job_context["success"])

        # Check output file
        output_file = self.test_dir + 'single_output/unmapped_by_salmon.fa'
        expected_output_file = self.test_dir + 'expected_single_output/unmapped_by_salmon.fa'
        self.assertTrue(identical_checksum(output_file, expected_output_file))


class TximportTestCase(TestCase):
    """Test salmon._tximport function, which launches tximport.R script."""

    def setUp(self):
        experiment = Experiment(accession_code='PRJNA408323')
        experiment.save()

        for id in ['07', '08', '09', '12', '13', '14']:
            sample = Sample(accession_code=('SRR60800' + id))
            sample.save()
            e_s = ExperimentSampleAssociation(experiment=experiment, sample=sample)
            e_s.save()

    @tag('salmon')
    def test_tximport_experiment(self):
        job_context = {
            'job_id': 456,
            'job': ProcessorJob(),
            'pipeline': Pipeline(name="Salmon"),
            'genes_to_transcripts_path': '/home/user/data_store/tximport_test/np_gene2txmap.txt'
        }

        experiment_dir = '/home/user/data_store/tximport_test/PRJNA408323'
        final_context = salmon._tximport(job_context, experiment_dir)

        expected_output_dir = '/home/user/data_store/tximport_test/expected_output'
        for filename in ['txi_out.RDS', 'gene_lengthScaledTPM.tsv']:
            output_path = experiment_dir + '/' + filename
            expected_output = expected_output_dir + '/' + filename
            self.assertTrue(identical_checksum(output_path, expected_output))

        # Check the individual files
        self.assertTrue(len(final_context['individual_files']), 6)
        for file in final_context['individual_files']:
            self.assertTrue(os.path.isfile(file.absolute_file_path))


class DetermineIndexLengthTestCase(TestCase):
    """Test salmon._determine_index_length function, which gets the salmon index length of a sample.
    For now, these tests only ensure that the output of the new faster salmon index function match
    that of the old one for the test data."""

    @tag('salmon')
    def test_salmon_determine_index_length_single_read(self):
        """Test that the right length is calculated when the sample has one read."""
        job, files = prepare_job()

        job_context = salmon._prepare_files({'original_files': [files[0]]})
        results = salmon._determine_index_length(job_context)

        self.assertEqual(results['index_length_raw'], 41)
        self.assertEqual(results['index_length'], 'short')

    @tag('salmon')
    def test_salmon_determine_index_length_double_read(self):
        """Test that the right length is calculated when the sample has two reads."""
        job, files = prepare_job()

        job_context = salmon._prepare_files({'original_files': files})
        results = salmon._determine_index_length(job_context)

        self.assertEqual(results['index_length_raw'], 41)
        self.assertEqual(results['index_length'], 'short')

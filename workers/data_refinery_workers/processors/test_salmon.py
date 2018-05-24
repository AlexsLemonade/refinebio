import hashlib
import os
import shutil
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
    Sample,
    ProcessorJobOriginalFileAssociation,
    OriginalFileSampleAssociation
)
from data_refinery_workers.processors import salmon, utils

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

    samp = Sample()
    samp.organism = homo_sapiens
    samp.save()

    computational_result = ComputationalResult()
    computational_result.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = homo_sapiens
    organism_index.result = computational_result
    organism_index.save()

    comp_file = ComputedFile()
    comp_file.absolute_file_path = "/home/user/data_store/processed/TEST/TRANSCRIPTOME_INDEX/Homo_sapiens_short.tar.gz"
    comp_file.result = computational_result
    comp_file.calculate_size()
    comp_file.calculate_sha1()
    comp_file.save()

    og_file = OriginalFile()
    og_file.source_filename = "ERR003000_1.fastq.gz"
    og_file.filename = "ERR003000_1.fastq.gz"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR003000_1.fastq.gz"
    og_file.save()

    og_file2 = OriginalFile()
    og_file2.source_filename = "ERR003000_2.fastq.gz"
    og_file2.filename = "ERR003000_2.fastq.gz"
    og_file2.absolute_file_path = "/home/user/data_store/raw/TEST/SALMON/ERR003000_2.fastq.gz"
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


def identical_checksum(file1, file2):
    """Confirm that the two files have identical checksum."""
    checksum_1 = hashlib.md5(open(file1, 'rb').read()).hexdigest()
    checksum_2 = hashlib.md5(open(file2, 'rb').read()).hexdigest()
    return checksum_1 == checksum_2


class SalmonTestCase(TestCase):

    @tag('salmon')
    def test_salmon(self):
        """ """
        # Ensure any computed files from previous tests are removed.
        try:
            os.remove("/home/user/data_store/raw/TEST/SALMON/processed/quant.sf")
        except FileNotFoundError:
            pass

        job, files = prepare_job()
        salmon.salmon(job.pk)
        job = ProcessorJob.objects.get(id=job.pk)
        self.assertTrue(job.success)

    def test_fastqc(self):

        job, og_files = prepare_job()
        win_context = {
            'job': job,
            'job_id': 789,
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

        job_context = {
            'job_id': 123,
            'input_file_path': self.test_dir + 'double_input/reads_1.fastq',
            'input_file_path_2': self.test_dir + 'double_input/reads_2.fastq',
            'output_directory': self.test_dir + 'double_output/'
        }
        job_context["job"] = ProcessorJob()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        samp = Sample()
        samp.organism = homo_sapiens
        samp.save()
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
        job_context = {
            'job_id': 456,
            'input_file_path': self.test_dir + 'single_input/single_read.fastq',
            'output_directory': self.test_dir + 'single_output/'
        }
        job_context["job"] = ProcessorJob()
        
        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        samp = Sample()
        samp.organism = homo_sapiens
        samp.save()
        job_context["sample"] = sample

        salmon._run_salmontools(job_context, False)

        # Confirm job status
        self.assertTrue(job_context["success"])

        # Check output file
        output_file = self.test_dir + 'single_output/unmapped_by_salmon.fa'
        expected_output_file = self.test_dir + 'expected_single_output/unmapped_by_salmon.fa'
        self.assertTrue(identical_checksum(output_file, expected_output_file))

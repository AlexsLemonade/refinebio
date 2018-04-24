import os
import shutil
from contextlib import closing
from django.test import TestCase
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
import pandas as pd

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

    return pj

class SalmonTestCase(TestCase):

    def test_salmon(self):
        """ """
        job = prepare_job()
        salmon.salmon(job.pk)

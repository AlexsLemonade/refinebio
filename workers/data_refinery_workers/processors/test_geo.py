import os
import shutil
from contextlib import closing
from django.test import TestCase
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_workers.processors import illumina, utils
import pandas as pd

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "ILLUMINA_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427%5Fnon%2Dnormalized%2Etxt.gz"
    og_file.filename = "GSE22427_non-normalized.txt"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/ILLUMINA/GSE22427_non-normalized.txt"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj

class IlluminaToPCLTestCase(TestCase):

    def test_illumina_to_pcl(self):
        """ """
        job = prepare_job()
        illumina.illumina_to_pcl(job.pk)

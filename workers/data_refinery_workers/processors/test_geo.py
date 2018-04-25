import os
from contextlib import closing
from django.test import TestCase
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    OriginalFileSampleAssociation
)
from data_refinery_workers.processors import illumina, agilent_twocolor, utils
import pandas as pd

def prepare_illumina_job():
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

    sample_names = ["LV-C&si-Control-1",
        "LV-C&si-Control-2",
        "LV-C&si-Control-3",
        "LV-C&si-EZH2-1",
        "LV-C&si-EZH2-2",
        "LV-C&si-EZH2-3",
        "LV-EZH2&si-EZH2-1",
        "LV-EZH2&si-EZH2-2",
        "LV-EZH2&si-EZH2-3",
        "LV-T350A&si-EZH2-1",
        "LV-T350A&si-EZH2-2",
        "LV-T350A&si-EZH2-3"
    ]

    for name in sample_names:
        sample = Sample()
        sample.accession_code = name
        sample.title = name
        sample.save()

        sa = SampleAnnotation()
        sa.sample = sample
        sa.data = {}
        sa.save()

        sample_assoc = OriginalFileSampleAssociation()
        sample_assoc.original_file = og_file
        sample_assoc.sample = sample
        sample_assoc.save()

    return pj

def prepare_agilent_twocolor_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AGILENT_TWOCOLOR_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE22900&format=file"
    og_file.filename = "GSM466597_95899_agilent.txt"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/AGILENT_TWOCOLOR/GSM466597_95899_agilent.txt"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj

class IlluminaToPCLTestCase(TestCase):

    def test_illumina_to_pcl(self):
        """ """
        job = prepare_illumina_job()
        illumina.illumina_to_pcl(job.pk)
        self.assertTrue(os.path.isfile('/home/user/data_store/raw/TEST/processed/GSE22427_non-normalized.PCL'))

class AgilentTwoColorTestCase(TestCase):

    def test_agilent_twocolor(self):
        """ """
        job = prepare_agilent_twocolor_job()
        agilent_twocolor.agilent_twocolor_to_pcl(job.pk)
        self.assertTrue(os.path.isfile('/home/user/data_store/raw/TEST/processed/GSM466597_95899_agilent.PCL'))

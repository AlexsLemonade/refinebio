import os
import shutil
from contextlib import closing
from django.test import TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_workers.processors import no_op, utils

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "NO_OP"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"
    og_file.filename = "GSM1426071_CD_colon_active_1.CEL"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/CEL/GSM1426071_CD_colon_active_1.CEL"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj

class NOOPTestCase(TestCase):

    @tag('no_op')
    def test_noop(self):
        """ """
        job = prepare_job()
        no_op.no_op_processor(job.pk)

    @tag('no_op')
    def test_convert_simple_pcl(self):
        """ """
        job = prepare_job()

        pj = ProcessorJob()
        pj.pipeline_applied = "NO_OP"
        pj.save()

        # ID_REF    VALUE
        og_file = OriginalFile()
        og_file.source_filename = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59079/E-GEOD-59079.processed.1.zip"
        og_file.filename = "GSM1426279_sample_table.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NOOP/GSM1426279_sample_table.txt"
        og_file.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = pj
        assoc1.save()

        no_op.no_op_processor(job.pk)

        # No header - ex 
        # AFFX-BioB-3_at  0.74218756
        og_file = OriginalFile()
        og_file.source_filename = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE10188&format=file"
        og_file.filename = "GSM269747-tbl-1.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NOOP/GSM269747-tbl-1.txt"
        og_file.save()
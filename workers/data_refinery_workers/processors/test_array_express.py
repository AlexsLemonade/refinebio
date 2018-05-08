import os
from django.test import TestCase, tag
from data_refinery_common.models import (
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_workers.processors import array_express

def prepare_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
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

class AffyToPCLTestCase(TestCase):

    @tag("affymetrix")
    def test_affy_to_pcl(self):
        """ """
        job = prepare_job()
        array_express.affy_to_pcl(job.pk)

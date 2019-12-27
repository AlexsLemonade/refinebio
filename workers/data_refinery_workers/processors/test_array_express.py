import os
import shutil
from django.test import TestCase, tag
from data_refinery_common.models import (
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    ComputationalResult,
    ComputedFile
)
from data_refinery_workers.processors import array_express


def prepare_ba_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"
    og_file.filename = "GSM1426071_CD_colon_active_1.CEL"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/CEL/GSM1426071_CD_colon_active_1.CEL"
    og_file.is_downloaded = True
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj

def prepare_non_ba_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM45nnn/GSM45588/suppl/GSM45588.CEL.gz"
    og_file.filename = "GSM45588.CEL"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/CEL/GSM45588.CEL"
    og_file.is_downloaded = True
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj


def prepare_huex_v1_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1364nnn/GSM1364667/suppl/GSM1364667_U_110208_7-02-10_S2.CEL.gz"
    og_file.filename = "GSM1364667_U_110208_7-02-10_S2.CEL"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/CEL/GSM1364667_U_110208_7-02-10_S2.CEL"
    og_file.is_downloaded = True
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
        job = prepare_ba_job()
        # Make sure that a previous test didn't leave a directory around.
        shutil.rmtree("/home/user/data_store/processor_job_" + str(job.id), ignore_errors=True)
        array_express.affy_to_pcl(job.pk)

        updated_job = ProcessorJob.objects.get(pk=job.pk)
        self.assertTrue(updated_job.success)
        self.assertEqual(len(ComputationalResult.objects.all()), 1)
        self.assertEqual(len(ComputedFile.objects.all()), 1)
        self.assertEqual(ComputedFile.objects.all()[0].filename, 'GSM1426071_CD_colon_active_1.PCL')

        os.remove(ComputedFile.objects.all()[0].absolute_file_path)

    @tag("affymetrix")
    def test_affy_to_pcl_no_brainarray(self):
        """ """
        job = prepare_non_ba_job()
        # Make sure that a previous test didn't leave a directory around.
        shutil.rmtree("/home/user/data_store/processor_job_" + str(job.id), ignore_errors=True)
        array_express.affy_to_pcl(job.pk)

        updated_job = ProcessorJob.objects.get(pk=job.pk)
        self.assertTrue(updated_job.success)
        self.assertEqual(len(ComputationalResult.objects.all()), 1)
        self.assertEqual(len(ComputedFile.objects.all()), 1)
        self.assertEqual(ComputedFile.objects.all()[0].filename, 'GSM45588.PCL')

        os.remove(ComputedFile.objects.all()[0].absolute_file_path)

    @tag("affymetrix", "huex")
    def test_affy_to_pcl_huex_v1(self):
        """ Special Case because there is no CDL for Huex V1 """
        job = prepare_huex_v1_job()
        shutil.rmtree("/home/user/data_store/processor_job_" + str(job.id), ignore_errors=True)
        array_express.affy_to_pcl(job.pk)

        updated_job = ProcessorJob.objects.get(pk=job.pk)
        self.assertTrue(updated_job.success)
        self.assertEqual(len(ComputationalResult.objects.all()), 1)
        self.assertEqual(len(ComputedFile.objects.all()), 1)
        self.assertEqual(ComputedFile.objects.all()[0].filename, 'GSM1364667_U_110208_7-02-10_S2.PCL')

        os.remove(ComputedFile.objects.all()[0].absolute_file_path)

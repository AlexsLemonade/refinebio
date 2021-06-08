import os
import shutil

from django.test import TestCase, tag

import pandas as pd
import scipy.stats

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
)
from data_refinery_workers.processors import array_express


def prepare_ba_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AFFY_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip"
    og_file.filename = "GSM1426071_CD_colon_active_1.CEL"
    og_file.absolute_file_path = (
        "/home/user/data_store/raw/TEST/CEL/GSM1426071_CD_colon_active_1.CEL"
    )
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
    og_file.source_filename = (
        "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM45nnn/GSM45588/suppl/GSM45588.CEL.gz"
    )
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
    og_file.absolute_file_path = (
        "/home/user/data_store/raw/TEST/CEL/GSM1364667_U_110208_7-02-10_S2.CEL"
    )
    og_file.is_downloaded = True
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj


class EnsgPkgMapTestCase(TestCase):
    @tag("affymetrix")
    def test_create_ensg_pkg_map(self):
        """Test _create_ensg_pkg_map using some arbitrary chips picked out of the list at
        http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/22.0.0/ensg.asp
        """

        ensg_pkg_map = array_express._create_ensg_pkg_map()
        self.assertIsNotNone(ensg_pkg_map)
        self.assertEquals(ensg_pkg_map["bovgene10st"], "bovgene10stbtensgprobe")
        self.assertEquals(ensg_pkg_map["hugene21st"], "hugene21sthsensgprobe")
        self.assertEquals(ensg_pkg_map["clariomdmouse"], "clariomdmousemmensgprobe")
        self.assertEquals(ensg_pkg_map["rattoxfx"], "rattoxfxrnensgprobe")
        self.assertEquals(ensg_pkg_map["yeast2"], "yeast2scensgprobe")


def assertMostlyAgrees(test_case: TestCase, expected_data: pd.Series, actual_data: pd.Series):
    """Checks to make sure that the expected data and the actual data are mostly
    the same, i.e. they have mostly the same genes and for the genes they have
    in common their spearman correlation is >0.99.

    We are only checking for approximate equality because output varies slightly
    between runs and between bioconductor versions."""

    # Make sure that the genes haven't changed too drastically between runs.
    # If this fails, it's probably not the end of the world but probably
    # something we should know about.
    test_case.assertGreater(
        len(set(expected_data.index) & set(actual_data.index)),
        0.95 * min(len(set(expected_data.index)), len(set(actual_data.index))),
    )

    expected_df = pd.DataFrame({"expected_values": expected_data})
    actual_df = pd.DataFrame({"actual_values": actual_data})

    # XXX: this doesn't work. The joining gets all messed up because some
    # indices have more than one associated value for some reason.
    (rho, _) = scipy.stats.spearmanr(expected_df.join(actual_df, how="inner"))
    test_case.assertGreater(rho, 0.99)


class AffyToPCLTestCase(TestCase):
    @tag("affymetrix")
    def test_affy_to_pcl(self):
        """ """
        job = prepare_ba_job()
        # Make sure that a previous test didn't leave a directory around.
        shutil.rmtree("/home/user/data_store/processor_job_" + str(job.id), ignore_errors=True)
        job_context = array_express.affy_to_pcl(job.pk)

        self.assertEqual(job_context["platform_accession_code"], "hugene10st")
        self.assertEqual(job_context["brainarray_package"], "hugene10sthsensgprobe")

        updated_job = ProcessorJob.objects.get(pk=job.pk)
        self.assertTrue(updated_job.success)
        self.assertEqual(len(ComputationalResult.objects.all()), 1)
        self.assertEqual(len(ComputedFile.objects.all()), 1)
        self.assertEqual(ComputedFile.objects.all()[0].filename, "GSM1426071_CD_colon_active_1.PCL")
        output_filename = ComputedFile.objects.all()[0].absolute_file_path

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/PCL/GSM1426071_CD_colon_active_1.PCL", sep="\t"
        )["GSM1426071_CD_colon_active_1.CEL"]
        actual_data = pd.read_csv(output_filename, sep="\t")["GSM1426071_CD_colon_active_1.CEL"]

        assertMostlyAgrees(self, expected_data, actual_data)

        os.remove(output_filename)

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
        self.assertEqual(ComputedFile.objects.all()[0].filename, "GSM45588.PCL")
        output_filename = ComputedFile.objects.all()[0].absolute_file_path

        expected_data = pd.read_csv("/home/user/data_store/TEST/PCL/GSM45588.PCL", sep="\t")[
            "GSM45588.CEL"
        ]
        actual_data = pd.read_csv(output_filename, sep="\t")["GSM45588.CEL"]

        assertMostlyAgrees(self, expected_data, actual_data)

        os.remove(output_filename)

    @tag("affymetrix", "huex")
    def test_affy_to_pcl_huex_v1(self):
        """Special Case because there is no CDL for Huex V1"""
        job = prepare_huex_v1_job()
        shutil.rmtree("/home/user/data_store/processor_job_" + str(job.id), ignore_errors=True)
        array_express.affy_to_pcl(job.pk)

        updated_job = ProcessorJob.objects.get(pk=job.pk)
        self.assertTrue(updated_job.success)
        self.assertEqual(len(ComputationalResult.objects.all()), 1)
        self.assertEqual(len(ComputedFile.objects.all()), 1)
        self.assertEqual(
            ComputedFile.objects.all()[0].filename, "GSM1364667_U_110208_7-02-10_S2.PCL"
        )
        output_filename = ComputedFile.objects.all()[0].absolute_file_path

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/PCL/GSM1364667_U_110208_7-02-10_S2.PCL", sep="\t"
        )["GSM1364667_U_110208_7-02-10_S2.CEL"]
        actual_data = pd.read_csv(output_filename, sep="\t")["GSM1364667_U_110208_7-02-10_S2.CEL"]

        assertMostlyAgrees(self, expected_data, actual_data)

        os.remove(ComputedFile.objects.all()[0].absolute_file_path)

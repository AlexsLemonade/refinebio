import math
import os
import shutil
from pathlib import Path

from django.test import TestCase, tag

import pandas as pd
import scipy.stats

from data_refinery_common.models import (
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
)
from data_refinery_workers.processors import no_op


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


class NOOPTestCase(TestCase):
    def setUp(self):
        """Cleanup work dirs from previous test runs that may have failed."""
        job_dirs = list(Path("/home/user/data_store").glob("**/processor_job_*"))
        for job_dir in job_dirs:
            shutil.rmtree(str(job_dir), ignore_errors=True)

    @tag("no_op")
    def test_convert_simple_pcl_with_header(self):
        """ """

        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        # ID_REF, VALUE
        og_file = OriginalFile()
        og_file.source_filename = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-51013/"
        og_file.filename = "GSM1234847_sample_table.txt"
        og_file.absolute_file_path = (
            "/home/user/data_store/raw/TEST/NO_OP/GSM1234847_sample_table.txt"
        )
        og_file.is_downloaded = True
        og_file.save()

        sample = Sample()
        sample.accession_code = "GSM1234847"
        sample.title = "GSM1234847"
        sample.platform_accession_code = "A-AFFY-38"
        sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.original_file = og_file
        assoc.sample = sample
        assoc.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        final_context = no_op.no_op_processor(job.pk)
        self.assertTrue(final_context["success"])
        self.assertTrue(os.path.exists(final_context["output_file_path"]))

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/NO_OP/EXPECTED/gene_converted_GSM1234847-tbl-1.txt",
            sep="\t",
            index_col=0,
        )["VALUE"]
        actual_data = pd.read_csv(final_context["output_file_path"], sep="\t", index_col=0)["VALUE"]

        assertMostlyAgrees(self, expected_data, actual_data)

    @tag("no_op")
    def test_convert_simple_pcl_with_no_header(self):
        """ """
        # No header - ex
        # AFFX-BioB-3_at  0.74218756
        og_file = OriginalFile()
        og_file.source_filename = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-51013/"
        og_file.filename = "GSM1234847_sample_table_headerless.txt"
        og_file.absolute_file_path = (
            "/home/user/data_store/raw/TEST/NO_OP/GSM1234847_sample_table_headerless.txt"
        )
        og_file.is_downloaded = True
        og_file.save()

        sample = Sample()
        sample.accession_code = "GSM1234847"
        sample.title = "GSM1234847"
        sample.platform_accession_code = "A-AFFY-38"
        sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.original_file = og_file
        assoc.sample = sample
        assoc.save()

        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        final_context = no_op.no_op_processor(job.pk)
        self.assertTrue(final_context["success"])
        self.assertTrue(os.path.exists(final_context["output_file_path"]))

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/NO_OP/EXPECTED/gene_converted_GSM1234847-tbl-1.txt",
            sep="\t",
            index_col=0,
        )["VALUE"]
        actual_data = pd.read_csv(final_context["output_file_path"], sep="\t", index_col=0)["VALUE"]

        assertMostlyAgrees(self, expected_data, actual_data)

    @tag("no_op")
    def test_convert_processed_illumina(self):
        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        # ex:
        # Reporter Identifier VALUE   Detection Pval
        # ILMN_1343291    14.943602   0
        # ILMN_1343295    13.528082   0
        og_file = OriginalFile()
        og_file.source_filename = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22433/"
        og_file.filename = "GSM557500-tbl-1.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM557500-tbl-1.txt"
        og_file.is_downloaded = True
        og_file.save()

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        homo_sapiens.save()

        sample = Sample()
        sample.accession_code = "GSM557500"
        sample.title = "GSM557500"
        sample.platform_accession_code = "A-MEXP-1171"
        sample.manufacturer = "ILLUMINA"
        sample.organism = homo_sapiens
        sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.original_file = og_file
        assoc.sample = sample
        assoc.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        # To:
        # ENSG00000156508 14.943602
        # ENSG00000111640 13.528082
        final_context = no_op.no_op_processor(job.pk)
        self.assertTrue(final_context["success"])
        self.assertTrue(os.path.exists(final_context["output_file_path"]))
        self.assertTrue(no_op.check_output_quality(final_context["output_file_path"]))

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/NO_OP/EXPECTED/gene_converted_GSM557500-tbl-1.txt",
            sep="\t",
            names=["", "VALUE"],
            index_col=0,
        )["VALUE"]
        actual_data = pd.read_csv(
            final_context["output_file_path"], sep="\t", names=["", "VALUE"], index_col=0
        )["VALUE"]

        assertMostlyAgrees(self, expected_data, actual_data)

    @tag("no_op")
    def test_convert_illumina_no_header(self):
        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        # ex:
        # ILMN_1885639    10.0000 0.7931
        # ILMN_2209417    10.0000 0.2029
        # ILMN_1765401    152.0873    0.0000
        og_file = OriginalFile()
        og_file.source_filename = (
            "https://github.com/AlexsLemonade/refinebio/files/2255178/GSM1089291-tbl-1.txt"
        )
        og_file.filename = "GSM1089291-tbl-1.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM1089291-tbl-1.txt"
        og_file.is_downloaded = True
        og_file.save()

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        homo_sapiens.save()

        sample = Sample()
        sample.accession_code = "GSM557500"
        sample.title = "GSM557500"
        sample.platform_accession_code = "A-MEXP-1171"
        sample.manufacturer = "ILLUMINA"
        sample.organism = homo_sapiens
        sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.original_file = og_file
        assoc.sample = sample
        assoc.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        # To:
        # ENSG00000105675 10
        # ENSG00000085721 152.0873
        # ENSG00000278494 152.0873
        final_context = no_op.no_op_processor(job.pk)
        self.assertTrue(final_context["success"])
        self.assertTrue(os.path.exists(final_context["output_file_path"]))

        expected_data = pd.read_csv(
            "/home/user/data_store/TEST/NO_OP/EXPECTED/gene_converted_GSM1089291-tbl-1.txt",
            sep="\t",
            names=["", "VALUE"],
            index_col=0,
        )["VALUE"]
        actual_data = pd.read_csv(
            final_context["output_file_path"], sep="\t", names=["", "VALUE"], index_col=0
        )["VALUE"]

        assertMostlyAgrees(self, expected_data, actual_data)

    @tag("no_op")
    def test_convert_illumina_bad_cols(self):
        """
        In future, this test may be deprecated. For now it just alerts that it needs attention.
        """
        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        # ex:
        # ILMN_1885639    10.0000 0.7931  11.0000 0.123
        # ILMN_2209417    10.0000 0.2029  11.1234 0.543
        # LMN_1765401    152.0873    0.0000  99.999  0.19
        og_file = OriginalFile()
        og_file.source_filename = (
            "https://github.com/AlexsLemonade/refinebio/files/2255178/GSM1089291-tbl-1-modified.txt"
        )
        og_file.filename = "GSM1089291-tbl-1-modified.txt"
        og_file.absolute_file_path = (
            "/home/user/data_store/raw/TEST/NO_OP/GSM1089291-tbl-1-modified.txt"
        )
        og_file.is_downloaded = True
        og_file.save()

        homo_sapiens = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        homo_sapiens.save()

        sample = Sample()
        sample.accession_code = "GSM557500"
        sample.title = "GSM557500"
        sample.platform_accession_code = "A-MEXP-1171"
        sample.manufacturer = "ILLUMINA"
        sample.organism = homo_sapiens
        sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.original_file = og_file
        assoc.sample = sample
        assoc.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = job
        assoc1.save()

        final_context = no_op.no_op_processor(job.pk)
        self.assertFalse(final_context["success"])
        self.assertTrue("Tell Rich!" in final_context["job"].failure_reason)

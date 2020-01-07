import os
import shutil
from pathlib import Path

from django.test import TestCase, tag

from data_refinery_common.models import (
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SurveyJob,
)
from data_refinery_workers.processors import no_op, utils


class NOOPTestCase(TestCase):
    def setUp(self):
        """Cleanup work dirs from previous test runs that may have failed."""
        job_dirs = list(Path("/home/user/data_store").glob("**/processor_job_*"))
        for job_dir in job_dirs:
            shutil.rmtree(str(job_dir), ignore_errors=True)

    @tag("no_op")
    def test_convert_simple_pcl(self):
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

        # No header - ex
        # AFFX-BioB-3_at  0.74218756
        og_file = OriginalFile()
        og_file.source_filename = (
            "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10188/miniml/GSE10188_family.xml.tgz"
        )
        og_file.filename = "GSM269747-tbl-1.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM269747-tbl-1.txt"
        og_file.is_downloaded = True
        og_file.save()

        sample = Sample()
        sample.accession_code = "GSM269747"
        sample.title = "GSM269747"
        sample.platform_accession_code = "GPL1319"
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
        self.assertEqual(os.path.getsize(final_context["output_file_path"]), 346535)

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
        og_file.filename = "GSM557500_sample_table.txt"
        og_file.absolute_file_path = (
            "/home/user/data_store/raw/TEST/NO_OP/GSM557500_sample_table.txt"
        )
        og_file.is_downloaded = True
        og_file.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

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
        self.assertEqual(os.path.getsize(final_context["output_file_path"]), 920374)
        self.assertTrue(no_op.check_output_quality(final_context["output_file_path"]))

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

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

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
        self.assertEqual(os.path.getsize(final_context["output_file_path"]), 786207)

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

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

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

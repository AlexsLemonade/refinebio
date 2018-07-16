import os
import shutil
from contextlib import closing
from django.test import TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    Sample,
    Organism,
    OriginalFileSampleAssociation,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_workers.processors import no_op, utils


def setUpModule():
    utils.createTestProcessors()


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
    def test_convert_simple_pcl(self):
        """ """

        job = ProcessorJob()
        job.pipeline_applied = "NO_OP"
        job.save()

        # ID_REF, VALUE
        og_file = OriginalFile()
        og_file.source_filename = "https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-51013/"
        og_file.filename = "GSM1234847_sample_table.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM1234847_sample_table.txt"
        og_file.save()

        sample = Sample()
        sample.accession_code = "GSM1234847"
        sample.title = "GSM1234847"
        sample.platform_accession_code = 'A-AFFY-38'
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
        og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10188/miniml/GSE10188_family.xml.tgz"
        og_file.filename = "GSM269747-tbl-1.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM269747-tbl-1.txt"
        og_file.save()

        sample = Sample()
        sample.accession_code = "GSM269747"
        sample.title = "GSM269747"
        sample.platform_accession_code = 'GPL1319'
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
        self.assertTrue(final_context['success'])

    @tag('no_op')
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
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/NO_OP/GSM557500_sample_table.txt"
        og_file.save()

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        sample = Sample()
        sample.accession_code = "GSM557500"
        sample.title = "GSM557500"
        sample.platform_accession_code = 'A-MEXP-1171'
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
        self.assertTrue(final_context['success'])
        self.assertTrue(os.path.exists(final_context['output_file_path']))

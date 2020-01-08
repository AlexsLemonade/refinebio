import os
import shutil
from unittest.mock import patch

from django.test import TestCase, tag

from data_refinery_common.models import (
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SurveyJob,
)
from data_refinery_workers.processors import transcriptome_index, utils


def prepare_job(length):

    pj = ProcessorJob()
    pj.pipeline_applied = "TRANSCRIPTOME_INDEX_" + length.upper()
    pj.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

    samp = Sample()
    samp.organism = homo_sapiens
    samp.accession_code = "derp" + length
    samp.save()

    og_file = OriginalFile()
    og_file.source_filename = "aegilops_tauschii_" + length + ".fa.gz"
    og_file.filename = "aegilops_tauschii_" + length + ".fa.gz"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/aegilops_tauschii_short.fa.gz"
    og_file.is_downloaded = True
    og_file.save()

    og_file2 = OriginalFile()
    og_file2.source_filename = "aegilops_tauschii_" + length + ".gtf.gz"
    og_file2.filename = "aegilops_tauschii_" + length + ".gtf.gz"
    og_file2.absolute_file_path = "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/aegilops_tauschii_short.gtf.gz"
    og_file2.is_downloaded = True
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

    assoc2 = ProcessorJobOriginalFileAssociation()
    assoc2.original_file = og_file2
    assoc2.processor_job = pj
    assoc2.save()

    return pj


class TXTestCase(TestCase):
    @tag("transcriptome")
    def test_assembly_information(self):
        og_file = OriginalFile()
        og_file.source_url = "ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz"
        og_file.source_filename = "aegilops_tauschii_short.fa.gz"

        og_file2 = OriginalFile()
        og_file2.source_url = "ftp://ftp.ensemblgenomes.org/pub/release-39/plants/gtf/aegilops_tauschii/Aegilops_tauschii.ASM34733v1.39.gtf.gz"
        og_file2.source_filename = "aegilops_tauschii_short.gtf.gz"

        job_context = {"original_files": [og_file, og_file2], "computed_files": []}
        job_context = transcriptome_index._extract_assembly_information(job_context)
        self.assertEqual("39", job_context["assembly_version"])
        self.assertEqual("ASM34733v1", job_context["assembly_name"])

    @tag("transcriptome")
    def test_tx(self):
        """ """
        # Make sure the work dirs don't exist cause this will fail the job.
        shutil.rmtree(
            "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/SHORT/processor_job_1",
            ignore_errors=True,
        )
        shutil.rmtree(
            "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/LONG/processor_job_2",
            ignore_errors=True,
        )

        job1 = prepare_job("short")
        job_context1 = transcriptome_index.build_transcriptome_index(job1.pk, length="short")
        job1 = ProcessorJob.objects.get(id=job1.pk)
        self.assertTrue(job1.success)
        self.assertEqual(job_context1["length"], "short")

        job2 = prepare_job("long")
        job_context2 = transcriptome_index.build_transcriptome_index(job2.pk, length="long")
        job2 = ProcessorJob.objects.get(id=job2.pk)
        self.assertTrue(job2.success)
        self.assertEqual(job_context2["length"], "long")

        self.assertNotEqual(job_context1["output_dir"], job_context2["output_dir"])

        self.assertTrue(os.path.exists(job_context1["computed_file"].get_synced_file_path()))
        self.assertTrue(os.path.exists(job_context2["computed_file"].get_synced_file_path()))
        self.assertNotEqual(
            job_context1["computed_file"].get_synced_file_path(),
            job_context2["computed_file"].get_synced_file_path(),
        )

        # This is the same logic as in `salmon._find_index`
        file = job_context1["computed_file"]
        unpacked = "/".join(file.get_synced_file_path().split("/")[:-1])
        self.assertTrue("SHORT" in unpacked)
        file2 = job_context2["computed_file"]
        unpacked2 = "/".join(file2.get_synced_file_path().split("/")[:-1])
        self.assertTrue("LONG" in unpacked2)

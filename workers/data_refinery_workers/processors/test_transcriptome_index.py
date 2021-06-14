import os
import shutil
from typing import List

from django.test import TestCase, tag

import pandas as pd

from data_refinery_common.job_lookup import ProcessorEnum
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


def prepare_original_files(length: str) -> List[OriginalFile]:
    og_file = OriginalFile()
    og_file.source_filename = "aegilops_tauschii_" + length + ".fa.gz"
    og_file.filename = "aegilops_tauschii_" + length + ".fa.gz"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/aegilops_tauschii_short.fa.gz"
    og_file.is_downloaded = True
    # We need to add the URL here so that _extract_assembly_information works properly
    og_file.source_url = "ftp://ftp.ensemblgenomes.org/pub/release-39/plants/fasta/aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz"
    og_file.save()

    og_file2 = OriginalFile()
    og_file2.source_filename = "aegilops_tauschii_" + length + ".gtf.gz"
    og_file2.filename = "aegilops_tauschii_" + length + ".gtf.gz"
    og_file2.absolute_file_path = "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/aegilops_tauschii_short.gtf.gz"
    og_file2.is_downloaded = True
    # We need to add the URL here so that _extract_assembly_information works properly
    og_file2.source_url = "ftp://ftp.ensemblgenomes.org/pub/release-39/plants/gtf/aegilops_tauschii/Aegilops_tauschii.ASM34733v1.39.gtf.gz"
    og_file2.save()

    return [og_file, og_file2]


def prepare_job(length):
    pj = ProcessorJob()
    pj.pipeline_applied = "TRANSCRIPTOME_INDEX_" + length.upper()
    pj.save()

    homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS", taxonomy_id=1001)

    samp = Sample()
    samp.organism = homo_sapiens
    samp.accession_code = "derp" + length
    samp.save()

    [og_file, og_file2] = prepare_original_files(length)

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
        job_context = {
            "computed_files": [],
            "length": "short",
            "cleanup": False,
        }
        job_context["original_files"] = prepare_original_files(job_context["length"])

        job_context = transcriptome_index._extract_assembly_information(job_context)
        self.assertEqual("39", job_context["assembly_version"])
        self.assertEqual("ASM34733v1", job_context["assembly_name"])
        self.assertEqual("EnsemblPlants", job_context["database_name"])

    @tag("transcriptome")
    def test_process_gtf(self):
        job_context = {
            "cleanup_gtf": False,
            "work_dir": "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII/SHORT/processor_job_process_gtf",
            "gtf_file_path": "/home/user/data_store/raw/TEST/TRANSCRIPTOME_INDEX/Homo_sapiens_testdata.gtf",
        }
        job_context["output_dir"] = job_context["work_dir"] + "/index"

        job_context = transcriptome_index._process_gtf(job_context)

        transcript_to_gene_ids = pd.read_csv(
            job_context["genes_to_transcripts_path"], sep="\t", index_col=1, names=["Gene"],
        )["Gene"]

        self.assertEqual(transcript_to_gene_ids.get("ENST00000480647", None), "ENSG00000041988")
        self.assertEqual(transcript_to_gene_ids.get("ENST00000486728", None), "ENSG00000186891")

        with open(job_context["gtf_file_path"], "r") as f:
            filtered_gtf = f.read()

        # ENST00000416931 is a pseudogene, so it should not appear in the filtered GTF
        self.assertNotIn('transcript_id "ENST00000416931"', filtered_gtf)
        # ENST00000480647 is not a pseudogene, so it should appear
        self.assertIn('transcript_id "ENST00000480647"', filtered_gtf)

    @tag("transcriptome")
    def test_get_organism_name_from_path(self):
        self.assertEqual(
            transcriptome_index.get_organism_name_from_path("/home/user/data_store/gallus_gallus"),
            "GALLUS_GALLUS",
        )

        self.assertEqual(
            transcriptome_index.get_organism_name_from_path(
                "/home/user/data_store/Homo_sapiens_foo_bar"
            ),
            "HOMO_SAPIENS",
        )

        self.assertEqual(
            transcriptome_index.get_organism_name_from_path(
                "/home/user/data_store/Candida_albicans_sc5314_gca_000784635"
            ),
            "CANDIDA_ALBICANS",
        )

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
        job_context1 = transcriptome_index.build_transcriptome_index(
            job1.pk, length="short", cleanup=False
        )
        job1 = ProcessorJob.objects.get(id=job1.pk)
        self.assertTrue(job1.success)
        self.assertEqual(job_context1["length"], "short")

        job2 = prepare_job("long")
        job_context2 = transcriptome_index.build_transcriptome_index(
            job2.pk, length="long", cleanup=False
        )
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


class RuntimeProcessorTest(TestCase):
    """Test the processor hosted inside the transcriptome index docker container."""

    @tag("transcriptome")
    def test_transcriptome(self):
        self.assertEqual(Processor.objects.count(), 0)  # No processor yet
        proc_key = "TX_INDEX"
        tx_processor = utils.find_processor(proc_key)
        self.assertEqual(Processor.objects.count(), 1)  # New processor created

        # Validate some information of the new processor
        self.assertEqual(tx_processor.name, ProcessorEnum[proc_key].value["name"])

        cmd_str = "salmon --version"
        cmd_output = utils.get_cmd_lines([cmd_str])[cmd_str]
        self.assertEqual(tx_processor.environment["cmd_line"][cmd_str], cmd_output)

        cmd_str = "rsem-calculate-expression --version"
        cmd_output = utils.get_cmd_lines([cmd_str])[cmd_str]
        self.assertEqual(tx_processor.environment["cmd_line"][cmd_str], cmd_output)

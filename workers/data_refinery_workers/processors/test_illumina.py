import os
import shutil
from typing import Dict

from django.test import TestCase, tag

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.models import (
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    Pipeline,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
)
from data_refinery_workers.processors import illumina, utils


def prepare_illumina_job(job_info: Dict) -> ProcessorJob:
    pj = ProcessorJob()
    pj.pipeline_applied = "ILLUMINA_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = job_info["source_filename"]
    og_file.filename = job_info["filename"]
    og_file.absolute_file_path = job_info["absolute_file_path"]
    og_file.is_downloaded = True
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    for s in job_info["samples"]:
        # For convenience, if you give a list of strings we'll just use the
        # strings as both titles and accessions.
        if type(s) == str:
            accession_code = s
            title = s
        else:
            accession_code, title = s

        sample = Sample()
        sample.accession_code = accession_code
        sample.title = title
        sample.organism = job_info["organism"]
        sample.save()

        sa = SampleAnnotation()
        sa.sample = sample
        sa.data = {"description": [title]}
        sa.is_ccdl = False
        sa.save()

        sample_assoc = OriginalFileSampleAssociation()
        sample_assoc.original_file = og_file
        sample_assoc.sample = sample
        sample_assoc.save()

    return pj


# Save this experiment separately (sans organism) because we need it for multiple tests
GSE22427 = {
    "source_filename": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427%5Fnon%2Dnormalized%2Etxt.gz",
    "filename": "GSE22427_non-normalized.txt",
    "absolute_file_path": "/home/user/data_store/raw/TEST/ILLUMINA/GSE22427_non-normalized.txt",
    "samples": [
        "LV-C&si-Control-1",
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
        "LV-T350A&si-EZH2-3",
    ],
}


class IlluminaToPCLTestCase(TestCase):
    @tag("illumina")
    def test_illumina_to_pcl(self):
        """Most basic Illumina to PCL test"""

        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        job = prepare_illumina_job({**GSE22427, "organism": organism})

        # Remove the title of one of the samples to make sure that we can still
        # find its detection column using the description given as an annotation
        sample = Sample.objects.get(title="LV-T350A&si-EZH2-3")
        sample.title = "ignoreme_for_description"
        sample.accession_code = "ignoreme_for_description"
        sample.save()

        final_context = illumina.illumina_to_pcl(job.pk)
        self.assertTrue(final_context["success"])

        for sample in final_context["samples"]:
            smashme = sample.get_most_recent_smashable_result_file()
            self.assertTrue(os.path.exists(smashme.absolute_file_path))
            os.remove(smashme.absolute_file_path)

        # Cleanup after the job since it won't since we aren't running in cloud.
        shutil.rmtree(final_context["work_dir"], ignore_errors=True)

    @tag("illumina")
    def test_bad_illumina_detection(self):
        """With the wrong species, this will fail the platform detection threshold."""

        organism = Organism(name="RATTUS_NORVEGICUS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        job = prepare_illumina_job({**GSE22427, "organism": organism})
        final_context = illumina.illumina_to_pcl(job.pk)
        self.assertTrue(final_context["abort"])

        # Cleanup after the job since it won't since we aren't running in cloud.
        shutil.rmtree(final_context["work_dir"], ignore_errors=True)

    @tag("illumina")
    def test_good_detection(self):
        """GSE54661 appears to be mislabled (illuminaHumanv4) on GEO. Shows our detector works."""

        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        pj = prepare_illumina_job(
            {
                "source_filename": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54661/suppl/GSE54661%5Fnon%5Fnormalized%2Etxt%2Egz",
                "filename": "GSE54661_non_normalized.txt",
                "absolute_file_path": "/home/user/data_store/raw/TEST/ILLUMINA/GSE54661_non_normalized.txt",
                "organism": organism,
                "samples": [("ABCD-1234", "CB CD34+ hypoxia"), ("ABCD-1235", "CB CD34+ normoxia")],
            }
        )

        final_context = illumina.illumina_to_pcl(pj.pk)
        self.assertTrue(final_context["success"])
        self.assertEqual(final_context["platform"], "illuminaHumanv3")

        for key in final_context["samples"][0].sampleannotation_set.all()[0].data.keys():
            self.assertTrue(
                key in ["detected_platform", "detection_percentage", "mapped_percentage"]
            )

        # Cleanup after the job since it won't since we aren't running in cloud.
        shutil.rmtree(final_context["work_dir"], ignore_errors=True)

    @tag("illumina")
    def test_illumina_latin1_input(self):
        """Test a latin1-encoded Illumina file.

        GSE106321 is encoded in latin1 and uses μ in the title of some
        columns, so preparing the file would cause a UnicodeParseError. Make
        sure that doesn't happen any more.
        """

        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        pj = prepare_illumina_job(
            {
                "source_filename": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106321/suppl/GSE106321_non-normalized.txt.gz",
                "filename": "GSE106321_non_normalized.txt",
                "absolute_file_path": "/home/user/data_store/raw/TEST/ILLUMINA/GSE106321_non-normalized.txt",
                "organism": organism,
                # NOTE: this isn't all the samples in the experiment. See below for why this doesn't matter
                "samples": [("GSM2835933", "A375 24h DMSO treatment")],
            }
        )

        # XXX: This experiment currently doesn't succeed processing for other
        # reasons, but the issue we are testing for only appears in
        # _prepare_files so we will run the pipeline up to that point.
        pipeline = Pipeline(name=PipelineEnum.ILLUMINA.value)
        final_context = utils.run_pipeline(
            {"job_id": pj.id, "pipeline": pipeline}, [utils.start_job, illumina._prepare_files,],
        )

        self.assertIsNone(final_context.get("success"))

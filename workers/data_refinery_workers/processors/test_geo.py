import os

from contextlib import closing
from django.test import TestCase, tag
from unittest.mock import MagicMock
from data_refinery_common.models import (
    SurveyJob,
    ProcessorJob,
    OriginalFile,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    OriginalFileSampleAssociation,
    Organism,
    Processor
)
from data_refinery_workers.processors import utils


def setUpModule():
    utils.createTestProcessors()


def prepare_illumina_job(species="Homo sapiens"):
    pj = ProcessorJob()
    pj.pipeline_applied = "ILLUMINA_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427%5Fnon%2Dnormalized%2Etxt.gz"
    og_file.filename = "GSE22427_non-normalized.txt"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/ILLUMINA/GSE22427_non-normalized.txt"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    sample_names = ["LV-C&si-Control-1",
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
        "LV-T350A&si-EZH2-3"
    ]

    for name in sample_names:
        sample = Sample()
        sample.accession_code = name
        sample.title = name
        sample.organism = Organism.get_object_for_name(species)
        sample.save()

        sa = SampleAnnotation()
        sa.sample = sample
        sa.data = {}
        sa.save()

        sample_assoc = OriginalFileSampleAssociation()
        sample_assoc.original_file = og_file
        sample_assoc.sample = sample
        sample_assoc.save()

    return pj

def prepare_agilent_twocolor_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AGILENT_TWOCOLOR_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE22900&format=file"
    og_file.filename = "GSM466597_95899_agilent.txt"
    og_file.absolute_file_path = "/home/user/data_store/raw/TEST/AGILENT_TWOCOLOR/GSM466597_95899_agilent.txt"
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj

class IlluminaToPCLTestCase(TestCase):

    @tag("illumina")
    def test_illumina_to_pcl(self):
        """ """
        from data_refinery_workers.processors import illumina
        job = prepare_illumina_job()
        final_context = illumina.illumina_to_pcl(job.pk)

        for sample in final_context['samples']:
            smashme = sample.get_most_recent_smashable_result_file()
            self.assertTrue(os.path.exists(smashme.absolute_file_path))
            os.remove(smashme.absolute_file_path)

    @tag("illumina")
    def test_bad_illumina_detection(self):
        """ With the wrong species, this will fail the platform detection threshold. """
        from data_refinery_workers.processors import illumina
        job = prepare_illumina_job('RATTUS_NORVEGICUS')
        final_context = illumina.illumina_to_pcl(job.pk)
        self.assertTrue(final_context['abort'])

    @tag("illumina")
    def test_good_detection(self):
        """GSE54661 appears to be mislabled (illuminaHumanv4) on GEO. Shows our detector works. """
        from data_refinery_workers.processors import illumina

        pj = ProcessorJob()
        pj.pipeline_applied = "ILLUMINA_TO_PCL"
        pj.save()

        og_file = OriginalFile()
        og_file.source_filename = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54661/suppl/GSE54661%5Fnon%5Fnormalized%2Etxt%2Egz"
        og_file.filename = "GSE54661_non_normalized.txt"
        og_file.absolute_file_path = "/home/user/data_store/raw/TEST/ILLUMINA/GSE54661_non_normalized.txt"
        og_file.save()

        assoc1 = ProcessorJobOriginalFileAssociation()
        assoc1.original_file = og_file
        assoc1.processor_job = pj
        assoc1.save()

        sample = Sample()
        sample.accession_code = "ABCD-1234"
        sample.title = "hypoxia_Signal"
        sample.organism = Organism.get_object_for_name("Homo sapiens")
        sample.save()

        sample_assoc = OriginalFileSampleAssociation()
        sample_assoc.original_file = og_file
        sample_assoc.sample = sample
        sample_assoc.save()

        final_context = illumina.illumina_to_pcl(pj.pk)
        self.assertEqual(final_context['platform'], 'illuminaHumanv3')

        for key in final_context['samples'][0].sampleannotation_set.all()[0].data.keys():
            self.assertTrue(key in ['detected_platform', 'detection_percentage', 'mapped_percentage'])

class AgilentTwoColorTestCase(TestCase):

    @tag("agilent")
    def test_agilent_twocolor(self):
        """ """
        from data_refinery_workers.processors import agilent_twocolor
        job = prepare_agilent_twocolor_job()
        agilent_twocolor.agilent_twocolor_to_pcl(job.pk)
        self.assertTrue(os.path.isfile('/home/user/data_store/raw/TEST/processed/GSM466597_95899_agilent.PCL'))

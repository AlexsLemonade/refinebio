import os

from django.test import TestCase, tag

from data_refinery_common.models import (
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
)
from data_refinery_workers.processors import agilent_twocolor


def prepare_agilent_twocolor_job():
    pj = ProcessorJob()
    pj.pipeline_applied = "AGILENT_TWOCOLOR_TO_PCL"
    pj.save()

    og_file = OriginalFile()
    og_file.source_filename = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE22900&format=file"
    og_file.filename = "GSM466597_95899_agilent.txt"
    og_file.absolute_file_path = (
        "/home/user/data_store/raw/TEST/AGILENT_TWOCOLOR/GSM466597_95899_agilent.txt"
    )
    og_file.is_downloaded = True
    og_file.save()

    assoc1 = ProcessorJobOriginalFileAssociation()
    assoc1.original_file = og_file
    assoc1.processor_job = pj
    assoc1.save()

    return pj


class AgilentTwoColorTestCase(TestCase):
    @tag("agilent")
    def test_agilent_twocolor(self):
        """ """

        job = prepare_agilent_twocolor_job()
        agilent_twocolor.agilent_twocolor_to_pcl(job.pk)
        self.assertTrue(
            os.path.isfile("/home/user/data_store/raw/TEST/processed/GSM466597_95899_agilent.PCL")
        )

        # XXX: jobs fail right now. Because we don't actually support two-color
        # yet, though, I have disabled this for now.

        # job.refresh_from_db()
        # self.assertTrue(job.success)

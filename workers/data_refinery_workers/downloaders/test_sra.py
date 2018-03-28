from urllib.error import URLError
from typing import List
from unittest.mock import patch, call
from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    OriginalFile,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob
)
from data_refinery_workers.downloaders import sra
from data_refinery_common.job_lookup import ProcessorPipeline


class DownloadSraTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        self.survey_job = survey_job

    def insert_objects(self):
        return 

    def test_download_file(self):
        dlj = DownloaderJob()
        dlj.accession_code = "ERR036"
        dlj.save()
        og = OriginalFile()
        og.source_filename = "ERR036000.fastq.gz"
        og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz"
        og.is_archive = True
        og.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og
        assoc.save()

        sra.download_sra(dlj.pk)

import os
from urllib.error import URLError
from typing import List
from unittest.mock import patch, call
from django.test import TestCase, tag
from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    OriginalFile,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob,
    Sample,
    OriginalFileSampleAssociation
)
from data_refinery_workers.downloaders import sra, utils
from data_refinery_common.job_lookup import ProcessorPipeline


class DownloadSraTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        self.survey_job = survey_job

    def insert_objects(self):
        return

    # @tag('downloaders')
    # @tag('downloaders_sra')
    # @patch('data_refinery_workers.downloaders.utils.send_job')
    # def test_download_file(self, mock_send_job):
    #     mock_send_job.return_value = None
        
    #     dlj = DownloaderJob()
    #     dlj.accession_code = "ERR036"
    #     dlj.save()

    #     og = OriginalFile()
    #     og.source_filename = "ERR036000.fastq.gz"
    #     og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz"
    #     og.is_archive = True
    #     og.save()

    #     sample = Sample()
    #     sample.accession_code = 'ERR036000'
    #     sample.save()

    #     assoc = OriginalFileSampleAssociation()
    #     assoc.sample = sample
    #     assoc.original_file = og
    #     assoc.save()

    #     assoc = DownloaderJobOriginalFileAssociation()
    #     assoc.downloader_job = dlj
    #     assoc.original_file = og
    #     assoc.save()

    #     success = sra.download_sra(dlj.pk)

    @tag('downloaders')
    @tag('downloaders_sra')
    @patch('data_refinery_workers.downloaders.utils.send_job')
    def test_download_file_ncbi(self, mock_send_job):
        mock_send_job.return_value = None
        
        dlj = DownloaderJob()
        dlj.accession_code = "DRR002116"
        dlj.save()
        og = OriginalFile()
        og.source_filename = "DRR002116.sra"
        og.source_url = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/DRR/DRR002/DRR002116/DRR002116.sra"
        og.is_archive = True
        og.save()
        sample = Sample()
        sample.accession_code = 'DRR002116'
        sample.save()
        assoc = OriginalFileSampleAssociation()
        assoc.sample = sample
        assoc.original_file = og
        assoc.save()
        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og
        assoc.save()
        result, downloaded_files = sra.download_sra(dlj.pk)
        utils.end_downloader_job(dlj, result)
        self.assertTrue(result)
        self.assertEqual(downloaded_files[0].sha1, 'd5374e7fe047d4f76b165c3f5148ab2df9d42cea')
        self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))
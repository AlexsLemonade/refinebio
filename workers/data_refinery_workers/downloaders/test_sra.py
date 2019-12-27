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
    OriginalFileSampleAssociation,
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

    @tag("downloaders")
    @tag("downloaders_sra")
    @patch("data_refinery_workers.downloaders.utils.send_job")
    def test_download_file_ncbi(self, mock_send_job):
        mock_send_job.return_value = None

        dlj = DownloaderJob()
        dlj.accession_code = "SRR9117853"
        dlj.save()
        og = OriginalFile()
        og.source_filename = "SRR9117853.sra"
        og.source_url = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR9117/SRR9117853/SRR9117853.sra"
        og.is_archive = True
        og.save()
        sample = Sample()
        sample.accession_code = "SRR9117853"
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
        self.assertEqual(
            downloaded_files[0].sha1, "e7ad484fe6f134ba7d1b2664e58cc15ae5a958cc"
        )
        self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))

    @tag("downloaders")
    @tag("downloaders_sra")
    @patch("data_refinery_workers.downloaders.utils.send_job")
    def test_download_file_swapper(self, mock_send_job):
        mock_send_job.return_value = None

        dlj = DownloaderJob()
        dlj.accession_code = "SRR9117853"
        dlj.save()
        og = OriginalFile()
        og.source_filename = "SRR9117853.sra"
        og.source_url = "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR9117/SRR9117853/SRR9117853.sra"
        og.is_archive = True
        og.save()
        sample = Sample()
        sample.accession_code = "SRR9117853"
        sample.save()
        assoc = OriginalFileSampleAssociation()
        assoc.sample = sample
        assoc.original_file = og
        assoc.save()
        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og
        assoc.save()
        result = sra._download_file(og.source_url, dlj, "/tmp/doomed", force_ftp=False)
        self.assertTrue(result)

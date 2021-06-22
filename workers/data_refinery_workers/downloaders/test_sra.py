import ftplib
import os
from ftplib import FTP
from unittest.mock import patch

from django.test import TestCase, tag

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    Sample,
    SurveyJob,
)
from data_refinery_workers.downloaders import sra, utils


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
    def test_download_file_ncbi(self):
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

        # If the FTP server is down or it times out, then we expect that the downloader job should have failed
        server_failure = False
        try:
            ftp_server = "ftp.sra.ebi.ac.uk"
            ftp = FTP(ftp_server)
            server_failure = ftp.login()[0:3] == "550"
        except (TimeoutError, ConnectionResetError):
            server_failure = True

        if server_failure:
            self.assertFalse(result)
        else:
            self.assertTrue(result)
            self.assertEqual(downloaded_files[0].sha1, "e7ad484fe6f134ba7d1b2664e58cc15ae5a958cc")
            self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))

    @tag("downloaders")
    @tag("downloaders_sra")
    @patch.object(ftplib.FTP, "login", side_effect=ftplib.error_perm)
    def test_ftp_server_down(self, mock_ftp):
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
        dlj.refresh_from_db()
        utils.end_downloader_job(dlj, result)
        self.assertFalse(result)
        self.assertEqual(dlj.failure_reason, "Failed to connect to ENA server.")

    @tag("downloaders")
    @tag("downloaders_sra")
    def test_download_file_swapper(self):
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

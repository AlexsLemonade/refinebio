import os
from unittest.mock import patch

from django.test import TestCase, tag

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
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

    @tag("downloaders")
    @tag("downloaders_sra")
    def test_download_file(self):
        dlj = DownloaderJob()
        dlj.accession_code = "ERR036"
        dlj.save()

        og = OriginalFile()
        og.source_filename = "ERR036000.fastq.gz"
        og.source_url = "ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz"
        og.is_archive = True
        og.save()

        sample = Sample()
        sample.accession_code = "ERR036000"
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

        self.assertTrue(result)
        self.assertEqual(downloaded_files[0].sha1, "1dfe5460a4101fe87feeffec0cb2e053f6695961")
        self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))

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

        self.assertTrue(result)
        self.assertEqual(downloaded_files[0].sha1, "e7ad484fe6f134ba7d1b2664e58cc15ae5a958cc")
        self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))

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
        result = sra._download_file(og, dlj, "/tmp/doomed")
        self.assertTrue(result)

    @tag("downloaders")
    @tag("downloaders_sra")
    def test_download_file_unmated_reads(self):
        dlj = DownloaderJob()
        dlj.accession_code = "SRR1603661"
        dlj.save()
        og_1 = OriginalFile()
        og_1.source_filename = "SRR1603661_1.fastq.gz"
        og_1.source_url = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR160/001/SRR1603661/SRR1603661_1.fastq.gz"
        og_1.expected_md5 = "502a9a482bfa5aa75865ccc0105ad13c"
        og_1.expected_size_in_bytes = 6751980628
        og_1.is_archive = True
        og_1.save()
        og_2 = OriginalFile()
        og_2.source_filename = "SRR1603661_2.fastq.gz"
        og_2.source_url = "ftp.sra.ebi.ac.uk/vol1/fastq/SRR160/001/SRR1603661/SRR1603661_2.fastq.gz"
        og_1.expected_md5 = "fffd24457418d255991f54ec82a39d57"
        og_1.expected_size_in_bytes = 6949912932
        og_2.is_archive = True
        og_2.save()
        sample = Sample()
        sample.accession_code = "SRR1603661"
        sample.save()
        assoc = OriginalFileSampleAssociation()
        assoc.sample = sample
        assoc.original_file = og_1
        assoc.save()
        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og_1
        assoc.save()
        assoc = OriginalFileSampleAssociation()
        assoc.sample = sample
        assoc.original_file = og_2
        assoc.save()
        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og_2
        assoc.save()
        result, downloaded_files = sra.download_sra(dlj.pk)
        utils.end_downloader_job(dlj, result)

        self.assertTrue(result)
        self.assertEqual(downloaded_files[0].sha1, "b5e227344deefad96dcc9ea209e5f6cc33be0b46")
        self.assertTrue(os.path.exists(downloaded_files[0].absolute_file_path))

import os
from django.test import TestCase, tag
from unittest.mock import patch

from data_refinery_common.models import (
    DownloaderJob,
    OriginalFile,
    DownloaderJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation
)
from data_refinery_workers.downloaders import geo, utils


class DownloadGeoTestCase(TestCase):
    def setUp(self):
        return

    @tag('downloaders')
    def test_download_and_extract_file(self):

        # Download function requires a DownloaderJob object,
        # can be blank for the simple case.
        dlj = DownloaderJob()
        dlj.save()

        # *_family.xml.tgz
        geo._download_file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10241/miniml/GSE10241_family.xml.tgz', '/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tgz', dlj)
        archive_file = geo.ArchivedFile('/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tgz')
        files = [file for file in archive_file.get_files()]

        # There should be 8 files in total in the directory, 2 downloaded and 6 extracted
        # `archive_file.get_files()` only returns the files that are extracted from the archives
        # instead of enumerating over all files.
        self.assertEqual(6, len(files))

        # GPL File
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GPL6102-tbl-1.txt'))

        # GSM Files
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSM258515-tbl-1.txt'))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSM258516-tbl-1.txt'))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSM258530-tbl-1.txt'))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSM258517-tbl-1.txt'))

        # Original family file
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSE10241_family.xml'))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tgz'))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tar'))

        # .txt.gz
        geo._download_file('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM254nnn/GSM254828/suppl/GSM254828.txt.gz', '/home/user/data_store/GSM254828/raw/GSM254828.txt.gz', dlj)
        archive_file = geo.ArchivedFile('/home/user/data_store/GSM254828/raw/GSM254828.txt.gz')
        files = [file for file in archive_file.get_files()]
        self.assertEqual(1, len(files))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSM254828/raw/GSM254828.txt'))

        geo._download_file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427_non-normalized.txt.gz", '/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt.gz',  dlj)
        archive_file = geo.ArchivedFile('/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt.gz')
        files = [file for file in archive_file.get_files()]
        self.assertEqual(1, len(files))

        self.assertTrue(os.path.isfile('/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt'))

    @tag('downloaders')
    def test_download_aspera_and_ftp(self):
        """ Tests the main 'download_geo' function. """

        dlj = DownloaderJob()
        dlj.accession_code = 'GSE22427'
        dlj.save()

        original_file = OriginalFile()
        original_file.source_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427_non-normalized.txt.gz"
        original_file.source_filename = "GSE22427_non-normalized.txt.gz"
        original_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.downloader_job = dlj
        assoc.save()

        sample = Sample()
        sample.accession_code = 'GSE22427'
        sample.save()

        sample_annotation = SampleAnnotation()
        sample_annotation.sample = sample
        sample_annotation.data = {'label_protocol_ch1': 'Agilent', 'label_protocol_ch2': 'Agilent'}
        sample_annotation.save()

        og_assoc = OriginalFileSampleAssociation()
        og_assoc.sample = sample
        og_assoc.original_file = original_file
        og_assoc.save()

        LOCAL_ROOT_DIR = "/home/user/data_store"
        os.makedirs(LOCAL_ROOT_DIR + '/' + sample.accession_code, exist_ok=True)
        dl_file_path = LOCAL_ROOT_DIR + '/' + sample.accession_code + '/' + original_file.source_url.split('/')[-1]

        # Aspera
        result = geo._download_file(original_file.source_url, file_path=dl_file_path, job=dlj, force_ftp=False)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(dl_file_path))
        os.remove(dl_file_path)

        # FTP
        result = geo._download_file(original_file.source_url, file_path=dl_file_path, job=dlj, force_ftp=True)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(dl_file_path))
        os.remove(dl_file_path)

        # Aspera, fail
        result = geo._download_file_aspera("https://rich.zone/cool_horse.jpg", target_file_path=dl_file_path, downloader_job=dlj, attempt=5)
        self.assertFalse(result)
        self.assertTrue(dlj.failure_reason != None)


    @tag('downloaders')
    @patch('data_refinery_workers.downloaders.utils.send_job')
    def test_download_geo(self, mock_send_task):
        """ Tests the main 'download_geo' function. """

        dlj = DownloaderJob()
        dlj.accession_code = 'GSE22427'
        dlj.save()

        original_file = OriginalFile()
        original_file.filename = "GSE22427_non-normalized.txt.gz"
        original_file.source_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427_non-normalized.txt.gz"
        original_file.source_filename = "GSE22427_non-normalized.txt.gz"
        original_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.downloader_job = dlj
        assoc.save()

        sample = Sample()
        sample.accession_code = 'GSE22427'
        sample.technology = "MICROARRAY"
        sample.manufacturer = "AGILENT"
        sample.has_raw = True
        # This is fake, but we don't currently support any agilent
        # platforms so we're using a platform that is supported.
        sample.platform_accession_code = "Illumina_RatRef-12_V1.0"
        sample.save()

        sample_annotation = SampleAnnotation()
        sample_annotation.sample = sample
        sample_annotation.data = {'label_protocol_ch1': 'Agilent', 'label_protocol_ch2': 'Agilent'}
        sample_annotation.save()

        og_assoc = OriginalFileSampleAssociation()
        og_assoc.sample = sample
        og_assoc.original_file = original_file
        og_assoc.save()

        download_result = geo.download_geo(dlj.id)

        file_assocs = OriginalFileSampleAssociation.objects.filter(sample=sample)
        self.assertEqual(file_assocs.count(), 2)

        for file_assoc in file_assocs:
            original_file = file_assoc.original_file
            if original_file.filename.endswith(".gz"):
                # We delete the archive after we extract from it
                self.assertFalse(original_file.is_downloaded)
            else:
                self.assertTrue(original_file.is_downloaded)

        # Make sure it worked
        self.assertTrue(download_result)
        self.assertTrue(dlj.failure_reason is None)
        self.assertTrue(len(ProcessorJob.objects.all()) > 0)
        self.assertEqual(ProcessorJob.objects.all()[0].pipeline_applied, "AGILENT_TWOCOLOR_TO_PCL")
        self.assertEqual(ProcessorJob.objects.all()[0].ram_amount, 2048)

    @tag('downloaders')
    def test_no_rnaseq(self):
        """Makes sure that no RNA-Seq data gets downloaded even if there's a job for it.
        """
        dlj = DownloaderJob()
        dlj.accession_code = 'GSE103217'
        dlj.save()

        original_file = OriginalFile()
        original_file.filename = "GSE103217_family.xml.tgz"
        original_file.source_url= "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103217/miniml/GSE103217_family.xml.tgz"
        original_file.source_filename = "GSE103217_family.xml.tgz"
        original_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.downloader_job = dlj
        assoc.save()

        sample = Sample()
        sample.accession_code = 'GSE103217'
        sample.technology = "RNA-SEQ"
        sample.manufacturer = "ILLUMINA"
        sample.platform_accession_code = "Illumina HiSeq 2500"
        sample.save()

        og_assoc = OriginalFileSampleAssociation()
        og_assoc.sample = sample
        og_assoc.original_file = original_file
        og_assoc.save()

        download_result = geo.download_geo(dlj.id)

        self.assertFalse(download_result)
        dlj.refresh_from_db()

        self.assertFalse(dlj.success)

        # It's not necessarily that we didn't extract any files, but
        # none that were usable so it looks like none.
        self.assertEqual(dlj.failure_reason, "Failed to extract any downloaded files.")

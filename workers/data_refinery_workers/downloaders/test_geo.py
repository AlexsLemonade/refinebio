import os

from django.test import TestCase
from data_refinery_common.models import (
    DownloaderJob,
    OriginalFile
)
from data_refinery_workers.downloaders import geo, utils


class DownloadGeoTestCase(TestCase):
    def setUp(self):
        return

    def test_download_and_extract_file(self):

        # Downlaod function requires a DownloaderJob object,
        # can be blank for the simple case.
        dlj = DownloaderJob()
        dlj.save()

        # *_family.xml.tgz
        geo._download_file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10241/miniml/GSE10241_family.xml.tgz', '/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tgz', dlj)
        files = geo._extract_tgz('/home/user/data_store/GSE10241/raw/GSE10241_family.xml.tgz', 'GSE10241')

        self.assertEqual(8, len(files))

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
        files = geo._extract_gz('/home/user/data_store/GSM254828/raw/GSM254828.txt.gz', 'GSM254828')
        self.assertEqual(1, len(files))
        self.assertTrue(os.path.isfile('/home/user/data_store/GSM254828/raw/GSM254828.txt'))

        geo._download_file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427%5Fnon%2Dnormalized%2Etxt%2Egz", '/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt.gz',  dlj)
        files = geo._extract_gz('/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt.gz', 'GSE22427')
        self.assertEqual(1, len(files))

        self.assertTrue(os.path.isfile('/home/user/data_store/GSE22427/raw/GSE22427_non-normalized.txt'))
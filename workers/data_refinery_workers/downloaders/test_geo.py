from django.test import TestCase
from data_refinery_common.models import (
    DownloaderJob,
    OriginalFile
)
from data_refinery_workers.downloaders import geo, utils


class DownloadGeoTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="GEO")
        survey_job.save()
        self.survey_job = survey_job

    def test_download_and_extract_file(self):
        dlj = DownloaderJob()
        dlj.save()

        # Miniml
        dlj = DownloaderJob()
        dlj.save()

        # .xml.tgz
        geo._download_file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10241/miniml/GSE10241_family.xml.tgz', 'GSE10241_family.xml.tgz', dlj)
        files = geo._extract_tgz('GSE10241_family.xml.tgz', 'GSE10241')

        # .txt.gz
        geo._download_file('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM254nnn/GSM254828/suppl/GSM254828.txt.gz', 'GSM254828.txt.gz', dlj)
        files = geo._extract_gz('GSM254828.txt.gz', 'GSM254828')
        geo._download_file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22427/suppl/GSE22427%5Fnon%2Dnormalized%2Etxt%2Egz", 'GSE22427.txt.gz',  dlj)


import copy
from typing import List
from unittest.mock import patch, call
from django.test import TestCase
from data_refinery_common.models import (
    SurveyJob,
    DownloaderJob,
    ProcessorJob
)
from data_refinery_workers.downloaders import array_express, utils
from data_refinery_common.job_lookup import ProcessorPipeline


class DownloadArrayExpressTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    @patch('data_refinery_workers.downloaders.utils.send_job')
    def test_download_and_extract_file(self, mock_send_job):
        dlj = DownloaderJob()
        dlj.save()
        array_express._download_file('ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip', 'dlme.zip', dlj)
        files = array_express._extract_files('dlme.zip', '123')


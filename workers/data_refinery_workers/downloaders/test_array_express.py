from unittest.mock import patch

from django.test import TestCase, tag

import vcr

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    Sample,
    SurveyJob,
)
from data_refinery_workers.downloaders import array_express


def file_caching_urlopen(download_url, *args, **kwargs):
    """Open a local file if it exists, otherwise download it and open it.

    Checks /home/user/data_store/test_download_files/registry.csv for
    `download_url`. If `download_url` is found in the registry open
    the file and return the file handler. If it is not found, download
    the file with urllib.request.urlopen(), add it to the registry,
    open the file, and return the file handler.
    """
    return open("/home/user/data_store/test_download_files/E-GEOD-59071.raw.3.zip", "rb")


class DownloadArrayExpressTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    @tag("downloaders")
    @patch("data_refinery_workers.downloaders.array_express.urllib.request.urlopen")
    @patch("data_refinery_workers.downloaders.utils.send_job")
    def test_download_and_extract_file(self, mock_send_job, mock_urlopen):
        mock_urlopen.side_effect = file_caching_urlopen
        dlj = DownloaderJob()
        dlj.save()
        array_express._download_file(
            "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",
            "dlme.zip",
            dlj,
        )
        files = array_express._extract_files("dlme.zip", "123", dlj)

    @tag("downloaders")
    @patch("data_refinery_workers.downloaders.utils.send_job")
    def test_download_multiple_zips(self, mock_send_job):
        """Tests that each sample gets one processor job no matter what.

        https://github.com/AlexsLemonade/refinebio/pull/351 deals with
        a bug where every file that was extracted to a directory got a
        processor job queued for it each time a downloader job ran
        which pointed to that directory. This test makes sure this bug
        stays squashed.

        It does so by running two downloader jobs for the same
        experiment which use two different zip files. Before this bug
        was squashed this would have resulted in the first sample
        getting a second processor job queued for it because the
        second downloader job would have found the file in the
        directory.
        """
        dlj1 = DownloaderJob()
        dlj1.accession_code = "E-MEXP-433"
        dlj1.save()

        original_file = OriginalFile()
        original_file.source_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MEXP/E-MEXP-433/E-MEXP-433.raw.1.zip"
        original_file.source_filename = "Waldhof_020604_R30_01-2753_U133A.CEL"
        original_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.downloader_job = dlj1
        assoc.save()

        sample = Sample()
        sample.accession_code = "E-MEXP-433-Waldhof_020604_R30_01-2753_U133A"
        sample.technology = "MICROARRAY"
        sample.manufacturer = "AFFYMETRIX"
        sample.has_raw = True
        # This is fake, but we don't currently support any agilent
        # platforms so we're using a platform that is supported.
        sample.platform_accession_code = "hgu133a"
        sample.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=sample, original_file=original_file
        )

        dlj2 = DownloaderJob()
        dlj2.accession_code = "E-MEXP-433"
        dlj2.save()

        original_file = OriginalFile()
        original_file.source_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MEXP/E-MEXP-433/E-MEXP-433.raw.2.zip"
        original_file.source_filename = "N08_U133A.CEL"
        original_file.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.downloader_job = dlj2
        assoc.save()

        sample = Sample()
        sample.accession_code = "E-MEXP-433-N08_U133A"
        sample.technology = "MICROARRAY"
        sample.manufacturer = "AFFYMETRIX"
        sample.has_raw = True
        # This is fake, but we don't currently support any agilent
        # platforms so we're using a platform that is supported.
        sample.platform_accession_code = "hgu133a"
        sample.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=sample, original_file=original_file
        )

        array_express.download_array_express(dlj1.id)
        array_express.download_array_express(dlj2.id)

        self.assertEqual(ProcessorJob.objects.all().count(), 2)

import csv
import os
import shutil
import urllib.request
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
from data_refinery_workers.downloaders import array_express

DOWNLOADED_FILES_DIR = "/home/user/data_store/test_download_files/"
DOWNLOADED_FILES_REGISTRY = DOWNLOADED_FILES_DIR + "test_file_registry.csv"
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256

original_urlopen = urllib.request.urlopen


def file_caching_urlopen(download_url, *args, **kwargs):
    """Open a local file if it exists, otherwise download it and open it.

    Checks /home/user/data_store/test_download_files/registry.csv for
    `download_url`. If `download_url` is found in the registry open
    the file and return the file handler. If it is not found, download
    the file with urllib.request.urlopen(), add it to the registry,
    open the file, and return the file handler.
    """
    registry = {}
    with open(DOWNLOADED_FILES_REGISTRY, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            registry[row["download_url"]] = row["file_path"]

    if download_url in registry:
        return open(registry[download_url], "rb")

    # We don't have it yet, we have to actually download it.
    file_path = DOWNLOADED_FILES_DIR + os.path.basename(download_url)
    with open(file_path, "wb") as target_file:
        with original_urlopen(download_url, timeout=60) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)

    with open(DOWNLOADED_FILES_REGISTRY, "a", newline="") as csvfile:
        fieldnames = ["download_url", "file_path"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writerow({"download_url": download_url, "file_path": file_path})

    return open(file_path, "rb")


class DownloadArrayExpressTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    @tag("downloaders")
    @patch("data_refinery_workers.downloaders.array_express.urllib.request.urlopen")
    def test_download_and_extract_file(self, mock_urlopen):
        mock_urlopen.side_effect = file_caching_urlopen
        dlj = DownloaderJob()
        dlj.save()
        array_express._download_file(
            "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEOD/E-GEOD-59071/E-GEOD-59071.raw.3.zip",
            "dlme.zip",
            dlj,
        )
        files = array_express._extract_files("dlme.zip", "123", dlj)

        # Test that all files were correctly extracted
        filenames = [file["filename"] for file in files]
        EXPECTED_FILES = [
            "GSM1426089_controle_colon_86.CEL",
            "GSM1426088_controle_colon_85.CEL",
            "GSM1426087_controle_colon_84.CEL",
            "GSM1426086_controle_colon_83.CEL",
            "GSM1426085_controle_colon_82.CEL",
            "GSM1426084_controle_colon_81.CEL",
            "GSM1426083_controle_colon_80.CEL",
            "GSM1426082_controle_colon_79.CEL",
            "GSM1426081_controle_colon_78.CEL",
            "GSM1426080_controle_colon_77.CEL",
            "GSM1426079_controle_colon_76.CEL",
            "GSM1426078_CD_colon_active_8.CEL",
            "GSM1426077_CD_colon_active_7.CEL",
            "GSM1426076_CD_colon_active_6.CEL",
            "GSM1426075_CD_colon_active_5.CEL",
            "GSM1426074_CD_colon_active_4.CEL",
            "GSM1426073_CD_colon_active_3.CEL",
            "GSM1426072_CD_colon_active_2.CEL",
            "GSM1426071_CD_colon_active_1.CEL",
        ]
        self.assertEqual(sorted(filenames), sorted(EXPECTED_FILES))

    @tag("downloaders")
    def test_download_multiple_zips(self):
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

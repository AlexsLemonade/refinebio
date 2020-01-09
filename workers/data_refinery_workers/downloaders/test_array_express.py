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


class DownloadArrayExpressTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="ARRAY_EXPRESS")
        survey_job.save()
        self.survey_job = survey_job

    @tag("downloaders")
    @patch("data_refinery_workers.downloaders.utils.send_job")
    def test_download_and_extract_file(self, mock_send_job):
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

    # @tag('downloaders')
    # def test_dharma(self):

    #     dlj1 = DownloaderJob()
    #     dlj1.accession_code = 'D1'
    #     dlj1.worker_id = get_instance_id()
    #     dlj1.start_time = datetime.datetime.now()
    #     dlj1.save()

    #     dlj2 = DownloaderJob()
    #     dlj2.accession_code = 'D2'
    #     dlj2.worker_id = get_instance_id()
    #     dlj2.start_time = datetime.datetime.now()
    #     dlj2.save()

    #     dlj3 = DownloaderJob()
    #     dlj3.accession_code = 'D3'
    #     dlj3.worker_id = get_instance_id()
    #     dlj3.save()

    #     original_file = OriginalFile()
    #     original_file.source_url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MEXP/E-MEXP-433/E-MEXP-433.raw.1.zip"
    #     original_file.source_filename = "Waldhof_020604_R30_01-2753_U133A.CEL"
    #     original_file.save()

    #     assoc = DownloaderJobOriginalFileAssociation()
    #     assoc.original_file = original_file
    #     assoc.downloader_job = dlj3
    #     assoc.save()

    #     sample = Sample()
    #     sample.accession_code = 'Blahblahblah'
    #     sample.technology = "MICROARRAY"
    #     sample.manufacturer = "AFFYMETRIX"
    #     sample.has_raw = True
    #     sample.platform_accession_code = "hgu133a"
    #     sample.save()

    #     OriginalFileSampleAssociation.objects.get_or_create(sample=sample, original_file=original_file)

    #     exited = False
    #     try:
    #         utils.start_job(dlj3.id, max_downloader_jobs_per_node=2, force_harakiri=True)
    #     except SystemExit as e:
    #         # This is supposed to happen!
    #         self.assertTrue(True)
    #         exited = True
    #     except Exception as e:
    #         # This isn't!
    #         self.assertTrue(False)
    #     self.assertTrue(exited)

    #     exited = False
    #     try:
    #         utils.start_job(dlj3.id, max_downloader_jobs_per_node=15, force_harakiri=True)
    #     except SystemExit as e:
    #         # This is not supposed to happen!
    #         self.assertTrue(False)
    #         exited = True
    #     except Exception as e:
    #         # This is!
    #         self.assertTrue(True)
    #     self.assertFalse(exited)

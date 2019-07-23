
from django.test import TransactionTestCase
from unittest.mock import Mock, patch, call
from django.test import tag

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.foreman.management.commands.create_missing_downloader_jobs import Command
from data_refinery_foreman.surveyor.geo import GeoSurveyor

class SurveyTestCase(TransactionTestCase):
    @tag('missing_jobs')
    def test_create_missing_jobs(self):
        """Tests that files which should have downloader jobs get them created."""

        # 1. create a sample with an original file and a downloader job
        original_file_with_downloader = OriginalFile()
        original_file_with_downloader.filename = "processed.CEL"
        original_file_with_downloader.source_filename = "processed.CEL"
        original_file_with_downloader.is_downloaded = True
        original_file_with_downloader.is_archive = False
        original_file_with_downloader.save()

        sample_with_downloader = Sample()
        sample_with_downloader.accession_code = "MA_doesnt_need_processor"
        sample_with_downloader.technology = "MICROARRAY"
        sample_with_downloader.source_database = "GEO"
        sample_with_downloader.platform_accession_code = "bovine"
        sample_with_downloader.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=sample_with_downloader,
            original_file=original_file_with_downloader
        )

        downloader_job = DownloaderJob()
        downloader_job.success = True
        downloader_job.worker_id = "worker_1"
        downloader_job.volume_index = "1"
        downloader_job.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=downloader_job,
            original_file=original_file_with_downloader
        )

        # 2. create a sample with an original file and no downloader job
        original_file = OriginalFile()
        original_file.filename = "tarball.gz"
        original_file.source_filename = "tarball.gz"
        original_file.is_downloaded = True
        original_file.is_archive = True
        original_file.save()

        sample_no_downloader = Sample()
        sample_no_downloader.accession_code = "sample_no_downloader"
        sample_no_downloader.technology = "MICROARRAY"
        sample_no_downloader.source_database = "GEO"
        sample_no_downloader.platform_accession_code = "bovine" # must be a supported platform
        sample_no_downloader.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=sample_no_downloader,
            original_file=original_file
        )

        # 3. Setup is done, actually run the command.
        command = Command()
        command.handle()

        ## Test that a missing downloader job was created.
        self.assertEqual(
            1,
            DownloaderJobOriginalFileAssociation.objects.filter(
                original_file=original_file
            ).count()
        )

        ## Test that a downloader job that wasn't missing wasn't created.
        ## Of course, we created one in test setup, so we're really
        ## checking that it's still only 1.
        self.assertEqual(
            1,
            DownloaderJobOriginalFileAssociation.objects.filter(
                original_file=original_file_with_downloader
            ).count()
        )

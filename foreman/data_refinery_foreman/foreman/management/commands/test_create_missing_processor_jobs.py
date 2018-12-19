from django.test import TransactionTestCase
from unittest.mock import Mock, patch, call

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.foreman.management.commands.create_missing_processor_jobs import Command
from data_refinery_foreman.surveyor.geo import GeoSurveyor

class SurveyTestCase(TransactionTestCase):
    # @patch('data_refinery_foreman.surveyor.external_source.message_queue.send_job')
    def test_create_missing_jobs(self):
        """Tests that files which should have processor jobs get them created.

        Specifically files that fall into this category are files that
        had successful downloader jobs but for some reason do not have
        processor jobs. It's not yet known why this is happening, but
        part of this management command is logging about them to get a
        grasp of how many there are.

        We want this test to cover both Microarray and RNA-Seq. We
        also need to test both that files which need processor jobs
        have them created, but also that files which don't need them
        don't get them created.

        Therefore we need 4 original files:
          * Microarray needing processor job.
          * Microarray not needing processor job.
          * RNA-Seq needing processor job.
          * Microarray not needing processor job.

        However Microarray can have files which shouldn't get
        processor jobs, so we're going to make one of those as
        well. Also Microarray jobs can download multiple files which
        get a processor job each, so we're going to make an additional
        Microarray file and associate it with the same downloader job
        so we can make sure two processor jobs are created based on
        that one downloader job.
        """

        ma_og_doesnt_need_processor = OriginalFile()
        ma_og_doesnt_need_processor.filename = "processed.CEL"
        ma_og_doesnt_need_processor.is_downloaded = True
        ma_og_doesnt_need_processor.is_archive = False
        ma_og_doesnt_need_processor.save()

        ma_sample_doesnt_need_processor = Sample()
        ma_sample_doesnt_need_processor.accession_code = "MA_doesnt_need_processor"
        ma_sample_doesnt_need_processor.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_doesnt_need_processor,
            original_file=ma_og_doesnt_need_processor
        )

        ma_dl_job_doesnt_need_processor = DownloaderJob()
        ma_dl_job_doesnt_need_processor.success = True
        ma_dl_job_doesnt_need_processor.worker_id = "worker_1"
        ma_dl_job_doesnt_need_processor.volume_index = "1"
        ma_dl_job_doesnt_need_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_doesnt_need_processor,
            original_file=ma_og_doesnt_need_processor
        )

        ma_processor_job = ProcessorJob()
        ma_processor_job.success = True
        ma_processor_job.worker_id = "worker_1"
        ma_dl_job_doesnt_need_processor.volume_index = "1"
        ma_processor_job.save()

        ProcessorJobOriginalFileAssociation.objects.get_or_create(
            processor_job=ma_processor_job,
            original_file=ma_og_doesnt_need_processor
        )

        ma_og_needs_processor_1 = OriginalFile()
        ma_og_needs_processor_1.filename = "something.CEL"
        ma_og_needs_processor_1.is_downloaded = True
        ma_og_needs_processor_1.is_archive = False
        ma_og_needs_processor_1.save()

        ma_og_needs_processor_2 = OriginalFile()
        ma_og_needs_processor_2.filename = "something_else.CEL"
        ma_og_needs_processor_2.is_downloaded = True
        ma_og_needs_processor_2.is_archive = False
        ma_og_needs_processor_2.save()

        ma_og_archive = OriginalFile()
        ma_og_archive.filename = "tarball.gz"
        ma_og_archive.is_downloaded = True
        ma_og_archive.is_archive = True
        ma_og_archive.save()

        ma_sample_needs_processor_1 = Sample()
        ma_sample_needs_processor_1.accession_code = "MA_needs_processor_1"
        ma_sample_needs_processor_1.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_1,
            original_file=ma_og_needs_processor_1
        )
        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_1,
            original_file=ma_og_archive
        )

        ma_sample_needs_processor_2 = Sample()
        ma_sample_needs_processor_2.accession_code = "MA_needs_processor_2"
        ma_sample_needs_processor_2.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_2,
            original_file=ma_og_needs_processor_2
        )
        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_2,
            original_file=ma_og_archive
        )

        ma_dl_job_needs_processor = DownloaderJob()
        ma_dl_job_needs_processor.success = True
        ma_dl_job_needs_processor.worker_id = "worker_1"
        ma_dl_job_doesnt_need_processor.volume_index = "1"
        ma_dl_job_needs_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor,
            original_file=ma_og_needs_processor_1
        )
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor,
            original_file=ma_og_needs_processor_2
        )
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor,
            original_file=ma_og_archive
        )

        # Setup is done, actually run the command.
        command = Command()
        command.handle()

        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_needs_processor_1
            ).count()
        )
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_needs_processor_2
            ).count()
        )
        self.assertEqual(
            0,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_archive
            ).count()
        )

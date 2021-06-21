from django.test import TransactionTestCase

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

        Therefore we need at least 4 original files:
          * Microarray needing processor job.
          * Microarray not needing processor job.
          * RNA-Seq needing processor job.
          * RNA-Seq not needing processor job.

        However Microarray can have files which shouldn't get
        processor jobs, so we're going to make one of those as
        well. Also Microarray jobs can download multiple files which
        get a processor job each, so we're going to make an additional
        Microarray file and associate it with the same downloader job
        so we can make sure two processor jobs are created based on
        that one downloader job.
        """
        # Microarray File/Samples/Jobs
        ma_og_doesnt_need_processor = OriginalFile()
        ma_og_doesnt_need_processor.filename = "processed.CEL"
        ma_og_doesnt_need_processor.source_filename = "processed.CEL"
        ma_og_doesnt_need_processor.is_downloaded = True
        ma_og_doesnt_need_processor.is_archive = False
        ma_og_doesnt_need_processor.save()

        ma_sample_doesnt_need_processor = Sample()
        ma_sample_doesnt_need_processor.accession_code = "MA_doesnt_need_processor"
        ma_sample_doesnt_need_processor.technology = "MICROARRAY"
        ma_sample_doesnt_need_processor.source_database = "ARRAY_EXPRESS"
        ma_sample_doesnt_need_processor.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_doesnt_need_processor, original_file=ma_og_doesnt_need_processor
        )

        ma_dl_job_doesnt_need_processor = DownloaderJob()
        ma_dl_job_doesnt_need_processor.success = True
        ma_dl_job_doesnt_need_processor.worker_id = "worker_1"
        ma_dl_job_doesnt_need_processor.volume_index = "1"
        ma_dl_job_doesnt_need_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_doesnt_need_processor,
            original_file=ma_og_doesnt_need_processor,
        )

        ma_processor_job = ProcessorJob()
        ma_processor_job.success = True
        ma_processor_job.worker_id = "worker_1"
        ma_processor_job.volume_index = "1"
        ma_processor_job.save()

        ProcessorJobOriginalFileAssociation.objects.get_or_create(
            processor_job=ma_processor_job, original_file=ma_og_doesnt_need_processor
        )

        ma_og_needs_processor_1 = OriginalFile()
        ma_og_needs_processor_1.filename = "something.CEL"
        ma_og_needs_processor_1.source_filename = "something.CEL"
        ma_og_needs_processor_1.is_downloaded = True
        ma_og_needs_processor_1.is_archive = False
        ma_og_needs_processor_1.save()

        ma_og_needs_processor_2 = OriginalFile()
        ma_og_needs_processor_2.filename = "something_else.CEL"
        ma_og_needs_processor_2.source_filename = "something_else.CEL"
        ma_og_needs_processor_2.is_downloaded = True
        ma_og_needs_processor_2.is_archive = False
        ma_og_needs_processor_2.save()

        ma_og_archive = OriginalFile()
        ma_og_archive.filename = "tarball.gz"
        ma_og_archive.source_filename = "tarball.gz"
        ma_og_archive.is_downloaded = True
        ma_og_archive.is_archive = True
        ma_og_archive.save()

        ma_sample_needs_processor_1 = Sample()
        ma_sample_needs_processor_1.accession_code = "MA_needs_processor_1"
        ma_sample_needs_processor_1.technology = "MICROARRAY"
        ma_sample_needs_processor_1.source_database = "ARRAY_EXPRESS"
        ma_sample_needs_processor_1.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_1, original_file=ma_og_needs_processor_1
        )
        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_1, original_file=ma_og_archive
        )

        ma_sample_needs_processor_2 = Sample()
        ma_sample_needs_processor_2.accession_code = "MA_needs_processor_2"
        ma_sample_needs_processor_2.technology = "MICROARRAY"
        ma_sample_needs_processor_2.source_database = "ARRAY_EXPRESS"
        ma_sample_needs_processor_2.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_2, original_file=ma_og_needs_processor_2
        )
        OriginalFileSampleAssociation.objects.get_or_create(
            sample=ma_sample_needs_processor_2, original_file=ma_og_archive
        )

        ma_dl_job_needs_processor = DownloaderJob()
        ma_dl_job_needs_processor.success = True
        ma_dl_job_needs_processor.worker_id = "worker_1"
        ma_dl_job_doesnt_need_processor.volume_index = "1"
        ma_dl_job_needs_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor, original_file=ma_og_needs_processor_1
        )
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor, original_file=ma_og_needs_processor_2
        )
        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=ma_dl_job_needs_processor, original_file=ma_og_archive
        )

        # RNA-Seq File/Samples/Jobs
        rna_og_doesnt_need_processor = OriginalFile()
        rna_og_doesnt_need_processor.filename = "processed.fastq"
        rna_og_doesnt_need_processor.source_filename = "processed.fastq"
        rna_og_doesnt_need_processor.is_downloaded = True
        rna_og_doesnt_need_processor.is_archive = False
        rna_og_doesnt_need_processor.save()

        rna_sample_doesnt_need_processor = Sample()
        rna_sample_doesnt_need_processor.accession_code = "RNA_doesnt_need_processor"
        rna_sample_doesnt_need_processor.technology = "RNA-SEQ"
        rna_sample_doesnt_need_processor.source_database = "SRA"
        rna_sample_doesnt_need_processor.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=rna_sample_doesnt_need_processor, original_file=rna_og_doesnt_need_processor
        )

        rna_dl_job_doesnt_need_processor = DownloaderJob()
        rna_dl_job_doesnt_need_processor.success = True
        rna_dl_job_doesnt_need_processor.worker_id = "worker_1"
        rna_dl_job_doesnt_need_processor.volume_index = "1"
        rna_dl_job_doesnt_need_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=rna_dl_job_doesnt_need_processor,
            original_file=rna_og_doesnt_need_processor,
        )

        rna_processor_job = ProcessorJob()
        # Failed ProcessorJobs will be retried, so they still count.
        rna_processor_job.success = False
        rna_processor_job.worker_id = "worker_1"
        rna_dl_job_doesnt_need_processor.volume_index = "1"
        rna_processor_job.save()

        ProcessorJobOriginalFileAssociation.objects.get_or_create(
            processor_job=rna_processor_job, original_file=rna_og_doesnt_need_processor
        )

        rna_og_needs_processor = OriginalFile()
        rna_og_needs_processor.filename = "something.fastq"
        rna_og_needs_processor.source_filename = "something.fastq"
        rna_og_needs_processor.is_downloaded = True
        rna_og_needs_processor.is_archive = False
        rna_og_needs_processor.save()

        rna_sample_needs_processor = Sample()
        rna_sample_needs_processor.accession_code = "RNA_needs_processor"
        rna_sample_needs_processor.technology = "RNA-SEQ"
        rna_sample_needs_processor.source_database = "SRA"
        rna_sample_needs_processor.save()

        OriginalFileSampleAssociation.objects.get_or_create(
            sample=rna_sample_needs_processor, original_file=rna_og_needs_processor
        )

        rna_dl_job_needs_processor = DownloaderJob()
        rna_dl_job_needs_processor.success = True
        rna_dl_job_needs_processor.worker_id = "worker_1"
        rna_dl_job_doesnt_need_processor.volume_index = "1"
        rna_dl_job_needs_processor.save()

        DownloaderJobOriginalFileAssociation.objects.get_or_create(
            downloader_job=rna_dl_job_needs_processor, original_file=rna_og_needs_processor
        )

        # Setup is done, actually run the command.
        command = Command()
        command.handle()

        # Test Microarray was handled correctly.
        ## Test that a missing processor job was created.
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_needs_processor_1
            ).count(),
        )
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_needs_processor_2
            ).count(),
        )
        self.assertEqual(
            0,
            ProcessorJobOriginalFileAssociation.objects.filter(original_file=ma_og_archive).count(),
        )

        ## Test that a processor job that wasn't missing wasn't created.
        ## Of course, we created one in test setup, so we're really
        ## checking that it's still only 1.
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=ma_og_doesnt_need_processor
            ).count(),
        )

        # Test Microarray was handled correctly.
        ## Test that the missing processor job was created.
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=rna_og_needs_processor
            ).count(),
        )

        ## Test that a processor job that wasn't missing wasn't created.
        ## Of course, we created one in test setup, so we're really
        ## checking that it's still only 1.
        self.assertEqual(
            1,
            ProcessorJobOriginalFileAssociation.objects.filter(
                original_file=rna_og_doesnt_need_processor
            ).count(),
        )

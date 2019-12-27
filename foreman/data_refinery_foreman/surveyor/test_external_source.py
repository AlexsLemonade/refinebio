from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    Sample,
    SurveyJob,
)


class SraSurveyorTestCase(TestCase):
    @patch("data_refinery_foreman.surveyor.external_source.message_queue.send_job")
    def test_queue_downloader_jobs_for_original_files(self, mock_send_task):
        """Make sure that queue_downloader_jobs queues all expected Downloader
        jobs for a given experiment.
        """
        # First, create an experiment with two samples associated with it
        # and create two original files for each of those samples.
        experiment_object = Experiment()
        experiment_object.accession_code = "Experiment1"
        experiment_object.save()

        sample_object_1 = Sample()
        sample_object_1.accession_code = "Sample1"
        sample_object_1.platform_accession_code = "Illumina Genome Analyzer"
        sample_object_1.platform_accession_name = "Illumina Genome Analyzer"
        sample_object_1.technology = "RNA-SEQ"
        sample_object_1.manufacturer = "ILLUMINA"
        sample_object_1.source_database = "SRA"
        sample_object_1.save()
        sample_object_2 = Sample()
        sample_object_2.accession_code = "Sample2"
        sample_object_2.platform_accession_code = "Illumina Genome Analyzer"
        sample_object_2.platform_accession_name = "Illumina Genome Analyzer"
        sample_object_2.technology = "RNA-SEQ"
        sample_object_2.manufacturer = "ILLUMINA"
        sample_object_2.source_database = "SRA"
        sample_object_2.save()

        association = ExperimentSampleAssociation()
        association.experiment = experiment_object
        association.sample = sample_object_1
        association.save()

        association = ExperimentSampleAssociation()
        association.experiment = experiment_object
        association.sample = sample_object_2
        association.save()

        sample_1_original_files = []
        sample_2_original_files = []

        original_file = OriginalFile()
        original_file.source_url = "first_url"
        original_file.source_filename = "first_filename"
        original_file.is_downloaded = False
        original_file.has_raw = True
        original_file.save()
        sample_1_original_files.append(original_file)

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file
        original_file_sample_association.sample = sample_object_1
        original_file_sample_association.save()

        original_file = OriginalFile()
        original_file.source_url = "second_url"
        original_file.source_filename = "second_filename"
        original_file.is_downloaded = False
        original_file.has_raw = True
        original_file.save()
        sample_2_original_files.append(original_file)

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file
        original_file_sample_association.sample = sample_object_1
        original_file_sample_association.save()

        original_file = OriginalFile()
        original_file.source_url = "third_url"
        original_file.source_filename = "third_filename"
        original_file.is_downloaded = False
        original_file.has_raw = True
        original_file.save()
        sample_2_original_files.append(original_file)

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file
        original_file_sample_association.sample = sample_object_2
        original_file_sample_association.save()

        original_file = OriginalFile()
        original_file.source_url = "fourth_url"
        original_file.source_filename = "fourth_filename"
        original_file.is_downloaded = False
        original_file.has_raw = True
        original_file.save()
        sample_2_original_files.append(original_file)

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file
        original_file_sample_association.sample = sample_object_2
        original_file_sample_association.save()

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        surveyor = SraSurveyor(survey_job)

        surveyor.queue_downloader_job_for_original_files(
            sample_1_original_files, experiment_object.accession_code
        )
        surveyor.queue_downloader_job_for_original_files(
            sample_2_original_files, experiment_object.accession_code
        )

        self.assertEqual(DownloaderJob.objects.all().count(), 2)

    def test_no_repeat_jobs(self):
        """Make sure that queue_downloader_jobs queues all expected Downloader
        jobs for a given experiment.
        """
        # First, create an experiment with two samples associated with it
        # and create two original files for each of those samples.
        experiment_object = Experiment()
        experiment_object.accession_code = "Experiment1"
        experiment_object.save()

        sample_object = Sample()
        sample_object.accession_code = "Sample1"
        sample_object.platform_accession_code = "Illumina Genome Analyzer"
        sample_object.platform_accession_name = "Illumina Genome Analyzer"
        sample_object.technology = "RNA-SEQ"
        sample_object.manufacturer = "ILLUMINA"
        sample_object.source_database = "SRA"
        sample_object.save()

        original_file_1 = OriginalFile()
        original_file_1.source_url = "first_url"
        original_file_1.source_filename = "first_filename"
        original_file_1.is_downloaded = False
        original_file_1.has_raw = True
        original_file_1.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file_1
        original_file_sample_association.sample = sample_object
        original_file_sample_association.save()

        original_file_2 = OriginalFile()
        original_file_2.source_url = "second_url"
        original_file_2.source_filename = "second_filename"
        original_file_2.is_downloaded = False
        original_file_2.has_raw = True
        original_file_2.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file_2
        original_file_sample_association.sample = sample_object
        original_file_sample_association.save()

        dlj = DownloaderJob()
        dlj.save()

        DownloaderJobOriginalFileAssociation(
            downloader_job=dlj, original_file=original_file_1
        ).save()

        DownloaderJobOriginalFileAssociation(
            downloader_job=dlj, original_file=original_file_2
        ).save()

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        surveyor = SraSurveyor(survey_job)

        surveyor.queue_downloader_job_for_original_files(
            [original_file_1, original_file_2], experiment_object.accession_code
        )

        # We made one DownloaderJob in this test, so
        # queue_downloader_job_for_original_files didn't have anything
        # to do, so there should still be only one:
        self.assertEqual(1, DownloaderJob.objects.all().count())

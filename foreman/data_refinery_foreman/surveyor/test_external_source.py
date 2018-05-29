from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_common.models import (
    SurveyJob,
    Experiment,
    ExperimentSampleAssociation,
    Sample,
    OriginalFile,
    OriginalFileSampleAssociation,
    DownloaderJob

)

class SraSurveyorTestCase(TestCase):
    @patch('data_refinery_foreman.surveyor.external_source.send_job')
    def test_queue_downloader_jobs(self, mock_send_task):
        """Make sure that queue_downloader_jobs queues all expected Downloader
        jobs for a given experiment.
        """
        # First, create an experiment with two samples associated with it
        # and create two original files for each of those samples.
        experiment_object = Experiment()
        experiment_object.save()

        sample_object_1 = Sample()
        sample_object_1.accession_code = "Sample1"
        sample_object_1.save()
        sample_object_2 = Sample()
        sample_object_2.accession_code = "Sample2"
        sample_object_2.save()

        association = ExperimentSampleAssociation()
        association.experiment = experiment_object
        association.sample = sample_object_1
        association.save()

        association = ExperimentSampleAssociation()
        association.experiment = experiment_object
        association.sample = sample_object_2
        association.save()

        original_file = OriginalFile()
        original_file.source_url = "first_url"
        original_file.source_filename = "first_filename"
        original_file.is_downloaded = False
        original_file.has_raw = True
        original_file.save()

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

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.original_file = original_file
        original_file_sample_association.sample = sample_object_2
        original_file_sample_association.save()

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        surveyor = SraSurveyor(survey_job)
        surveyor.queue_downloader_jobs(experiment_object)

        self.assertEqual(DownloaderJob.objects.all().count(), 4)

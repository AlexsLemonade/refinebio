from unittest import mock

from django.core.management import call_command
from django.test import TestCase

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)


class TestTrackDataRequest(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        experiment = Experiment.objects.create(
            accession_code="GSE12417", technology="MICROARRAY", source_database="GEO"
        )

        cls.sample1 = Sample.objects.create(
            accession_code="GSM311750",
            source_database="GEO",
            technology="MICROARRAY",
            platform_accession_code="hgu133a",
        )
        cls.sample2 = Sample.objects.create(
            accession_code="GSM311751",
            source_database="GEO",
            technology="MICROARRAY",
            platform_accession_code="hgu133plus2",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=cls.sample1)
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=cls.sample2)

    def call_command(self, experiments=(), **kwargs):
        with mock.patch(
            f"data_refinery_foreman.foreman.management.commands.track_data_request.Command.get_accessions",
            return_value=experiments,
        ):
            call_command("track_data_request", (), **{"source_url": "https://some.url"})

    @mock.patch("builtins.print")
    def test_no_experiments(self, output):
        self.call_command(experiments=())
        output.assert_called_once_with("No experiments found")

    @mock.patch("builtins.print")
    def test_no_experiments_available(self, output):
        self.call_command(experiments=("GSE12417",))
        output.assert_called_once_with(
            "\n".join(
                (
                    "Experiments attempted: 1",
                    "GEO: 1",
                    "",
                    "Samples attempted: 2",
                    "GEO: 2",
                    "",
                    "Experiments available: 0",
                    "",
                    "Samples available: 0",
                    "",
                    "Total jobs: 0",
                )
            )
        )

    @mock.patch("builtins.print")
    def test_available(self, output):
        self.sample1.is_processed = True
        self.sample1.save()
        self.sample2.is_processed = True
        self.sample2.save()

        self.call_command(experiments=("GSE12417",))
        output.assert_called_once_with(
            "\n".join(
                (
                    "Experiments attempted: 1",
                    "GEO: 1",
                    "",
                    "Samples attempted: 2",
                    "GEO: 2",
                    "",
                    "Experiments available: 1",
                    "GEO: 1",
                    "",
                    "Samples available: 2",
                    "GEO: 2",
                    "",
                    "Total jobs: 0",
                )
            )
        )

    @mock.patch("builtins.print")
    def test_multiple(self, output):
        self.sample1.is_processed = True
        self.sample1.save()
        self.sample2.is_processed = True
        self.sample2.save()

        experiment = Experiment.objects.create(
            accession_code="GSE12418",
            technology="MICROARRAY",
            source_database="ARRAY_EXPRESS",
        )
        sample3 = Sample.objects.create(
            accession_code="GSM311760",
            is_processed=True,
            technology="MICROARRAY",
        )
        sample4 = Sample.objects.create(
            accession_code="GSM311761",
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample3)
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample4)

        survey_job = SurveyJob(source_type="GEO")
        survey_job.save()

        # Jobs.
        SurveyJobKeyValue(
            survey_job=survey_job,
            key="experiment_accession_code",
            value=experiment.accession_code,
        ).save()

        original_file = OriginalFile()
        original_file.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.sample = sample3
        original_file_sample_association.original_file = original_file
        original_file_sample_association.save()

        downloader_job = DownloaderJob()
        downloader_job.save()

        download_association = DownloaderJobOriginalFileAssociation()
        download_association.original_file = original_file
        download_association.downloader_job = downloader_job
        download_association.save()

        processor_job = ProcessorJob(downloader_job=downloader_job)
        processor_job.save()

        processor_association = ProcessorJobOriginalFileAssociation()
        processor_association.original_file = original_file
        processor_association.processor_job = processor_job
        processor_association.save()

        self.call_command(experiments=("GSE12417", "GSE12418"))
        output.assert_called_once_with(
            "\n".join(
                (
                    "Experiments attempted: 2",
                    "ARRAY_EXPRESS: 1, GEO: 1",
                    "",
                    "Samples attempted: 4",
                    "ARRAY_EXPRESS: 2, GEO: 2",
                    "",
                    "Experiments available: 2",
                    "ARRAY_EXPRESS: 1, GEO: 1",
                    "",
                    "Samples available: 3",
                    "ARRAY_EXPRESS: 1, GEO: 2",
                    "",
                    "Total jobs: 3",
                )
            )
        )

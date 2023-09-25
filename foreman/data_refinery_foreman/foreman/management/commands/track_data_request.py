import tempfile
from collections import Counter

from django.core.management.base import BaseCommand

import boto3
import botocore

from data_refinery_common.models import DownloaderJob, Experiment, ProcessorJob, SurveyJob
from data_refinery_common.utils import parse_s3_url


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--file",
            type=str,
            help="A file listing accession codes. s3:// URLs only are accepted.",
        )

    def handle(self, *args, **options):
        """Creates a report based on the processing results."""
        experiments_attempted = 0
        experiments_attempted_by_source = Counter()
        experiments_available = 0
        experiments_available_by_source = Counter()
        jobs_created = 0
        jobs_created_by_source = Counter()
        samples_attempted = 0
        samples_attempted_by_source = Counter()
        samples_available = 0
        samples_available_by_source = Counter()

        accessions = self.get_accessions(options["file"])
        for accession in accessions or ():
            experiment = Experiment.objects.filter(accession_code=accession).first()
            if not experiment:
                continue

            source = experiment.source_database

            # Experiments attempted, total and breakdown by source.
            experiments_attempted += 1
            experiments_attempted_by_source.update({source: 1})

            # Samples attempted, total and breakdown by source.
            samples_count = experiment.samples.count()
            samples_attempted += samples_count
            samples_attempted_by_source.update({source: samples_count})

            # Experiments available, total and breakdown by source.
            processed_samples_count = experiment.samples.filter(is_processed=True).count()
            if processed_samples_count > 0:
                experiments_available += 1
                experiments_available_by_source.update({source: 1})

                # Samples available, total and breakdown by source.
                samples_available += processed_samples_count
                samples_available_by_source.update({source: processed_samples_count})

            sample_accessions = experiment.samples.values_list("accession_code", flat=True)
            # Total number of jobs created, breakdown by source.
            downloader_jobs_created = self.get_downloader_jobs(sample_accessions).count()
            processor_jobs_created = self.get_processor_jobs(sample_accessions).count()
            survey_jobs_created = self.get_surveyor_jobs(experiment).count()
            total_jobs_created = (
                downloader_jobs_created + processor_jobs_created + survey_jobs_created
            )

            jobs_created += total_jobs_created
            jobs_created_by_source.update({source: total_jobs_created})

            # TODO(arkid15r): Calculate total run time as a difference between first created job
            # and the last finished job. Also indicate if there is something that still needs to
            # be processed.

        output = []
        if experiments_attempted:
            output += [
                f"Experiments attempted: {experiments_attempted}",
                f"{self.get_distribution_by_source(experiments_attempted_by_source)}",
                "",
            ]

            output.append(f"Samples attempted: {samples_attempted}")
            if samples_attempted:
                output.append(f"{self.get_distribution_by_source(samples_attempted_by_source)}")
            output.append("")

            output.append(f"Experiments available: {experiments_available}")
            if experiments_available:
                output.append(f"{self.get_distribution_by_source(experiments_available_by_source)}")
            output.append("")

            output.append(f"Samples available: {samples_available}")
            if samples_available:
                output.append(f"{self.get_distribution_by_source(samples_available_by_source)}")
            output.append("")

            output.append(f"Total jobs: {total_jobs_created}")
        elif accessions is not None:
            output.append("No experiments found")

        if output:
            print("\n".join(output))

    @staticmethod
    def get_accessions(s3_url):
        """Gets source experiment accessions."""
        if not s3_url.startswith("s3://"):
            print("Please provide a valid S3 URL")
            return None

        with tempfile.TemporaryFile() as tmp_file:
            bucket, key = parse_s3_url(s3_url)
            try:
                boto3.resource("s3").meta.client.download_fileobj(bucket, key, tmp_file)
            except botocore.exceptions.ClientError as e:
                if e.response["Error"]["Code"] == "404":
                    print("The S3 file does not exist.")
                    return None

            return (line.strip() for line in tmp_file.readlines() if line.strip())

    @staticmethod
    def get_downloader_jobs(sample_accessions):
        """Returns downloader jobs for sample accessions."""
        return DownloaderJob.objects.filter(
            original_files__samples__accession_code__in=sample_accessions
        ).distinct()

    @staticmethod
    def get_processor_jobs(sample_accessions):
        """Returns processor jobs for sample accessions."""
        return ProcessorJob.objects.filter(
            original_files__samples__accession_code__in=sample_accessions
        ).distinct()

    @staticmethod
    def get_surveyor_jobs(experiment):
        """Returns surveyor jobs for an experiment."""
        return SurveyJob.objects.filter(
            surveyjobkeyvalue__key="experiment_accession_code",
            surveyjobkeyvalue__value=experiment.accession_code,
        )

    @staticmethod
    def get_distribution_by_source(stats):
        """Returns a source based stats in a consistent manner."""
        return ", ".join((f"{source}: {stats[source]}" for source in sorted(stats)))

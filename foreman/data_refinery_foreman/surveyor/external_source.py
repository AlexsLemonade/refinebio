import abc
import os

from django.db import transaction
from retrying import retry
from typing import List, Dict

from data_refinery_common import message_queue, job_lookup, logging
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentSampleAssociation,
    OriginalFile,
    Sample,
    SurveyJob,
)


logger = logging.get_and_configure_logger(__name__)


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    @abc.abstractmethod
    def discover_experiments_and_samples(self):
        """Abstract method to survey a source."""
        return

    def queue_downloader_jobs(self, experiment: Experiment, samples: List[Sample]):
        """This enqueues DownloaderJobs on a per-file basis.

        There is a complementary function below for enqueueing multi-file
        DownloaderJobs.
        """
        files_to_download = []
        for sample in samples:
            files_for_sample = OriginalFile.objects.filter(sample=sample, is_downloaded=False)
            for og_file in files_for_sample:
                files_to_download.append(og_file)

        download_urls_with_jobs = {}
        for original_file in files_to_download:

            # We don't need to create multiple downloaders for the same file.
            # However, we do want to associate original_files with the
            # DownloaderJobs that will download them.
            if original_file.source_url in download_urls_with_jobs.keys():
                DownloaderJobOriginalFileAssociation.objects.get_or_create(
                    downloader_job=download_urls_with_jobs[original_file.source_url],
                    original_file=original_file,
                )
                continue

            # There is already a downloader job associated with this file.
            old_assocs_count = DownloaderJobOriginalFileAssociation.objects.filter(
                original_file__source_url=original_file.source_url
            ).count()
            if old_assocs_count > 0:
                logger.debug(
                    "We found an existing DownloaderJob for this file/url.",
                    original_file_id=original_file.id,
                )
                continue

            sample_object = original_file.samples.first()
            downloader_task = job_lookup.determine_downloader_task(sample_object)

            if downloader_task == job_lookup.Downloaders.NONE:
                logger.info(
                    "No valid downloader task found for sample.",
                    sample=sample_object.id,
                    original_file=original_file.id,
                )
            else:
                downloader_job = DownloaderJob()
                downloader_job.downloader_task = downloader_task.value
                downloader_job.accession_code = experiment.accession_code
                downloader_job.save()

                DownloaderJobOriginalFileAssociation.objects.get_or_create(
                    downloader_job=downloader_job, original_file=original_file
                )

                download_urls_with_jobs[original_file.source_url] = downloader_job

                try:
                    logger.info(
                        "Queuing downloader job for URL: " + original_file.source_url,
                        survey_job=self.survey_job.id,
                        downloader_job=downloader_job.id,
                    )
                    message_queue.send_job(downloader_task, downloader_job)
                except Exception as e:
                    # If the task doesn't get sent we don't want the
                    # downloader_job to be left floating
                    logger.exception(
                        "Failed to enqueue downloader job for URL: " + original_file.source_url,
                        survey_job=self.survey_job.id,
                        downloader_job=downloader_job.id,
                    )
                    downloader_job.success = False
                    downloader_job.failure_reason = str(e)
                    downloader_job.save()

    def queue_downloader_job_for_original_files(
        self,
        original_files: List[OriginalFile],
        experiment_accession_code: str = None,
        is_transcriptome: bool = False,
    ):
        """Creates a single DownloaderJob with multiple files to download.
        """
        # Transcriptome is a special case because there's no sample_object.
        # It's alright to re-process transcriptome indices.
        if is_transcriptome:
            downloader_task = job_lookup.Downloaders.TRANSCRIPTOME_INDEX
        else:
            source_urls = [original_file.source_url for original_file in original_files]
            # There is already a downloader job associated with this file.
            old_assocs_count = DownloaderJobOriginalFileAssociation.objects.filter(
                original_file__source_url__in=source_urls
            ).count()
            if old_assocs_count > 0:
                logger.debug(
                    "We found an existing DownloaderJob for these urls.", source_urls=source_urls
                )
                return False

            sample_object = original_files[0].samples.first()
            downloader_task = job_lookup.determine_downloader_task(sample_object)

        if downloader_task == job_lookup.Downloaders.NONE:
            logger.info(
                "No valid downloader task found for sample.",
                sample=sample_object.id,
                original_file=original_files[0].id,
            )
        else:
            downloader_job = DownloaderJob()
            downloader_job.downloader_task = downloader_task.value
            downloader_job.accession_code = experiment_accession_code
            downloader_job.save()

            downloaded_urls = []
            for original_file in original_files:
                DownloaderJobOriginalFileAssociation.objects.get_or_create(
                    downloader_job=downloader_job, original_file=original_file
                )

                downloaded_urls.append(original_file.source_url)

            try:
                logger.info(
                    "Queuing downloader job.",
                    survey_job=self.survey_job.id,
                    downloader_job=downloader_job.id,
                    downloaded_urls=downloaded_urls,
                )
                message_queue.send_job(downloader_task, downloader_job)
            except Exception as e:
                # If the task doesn't get sent we don't want the
                # downloader_job to be left floating
                logger.exception(
                    "Failed to enqueue downloader job.",
                    survey_job=self.survey_job.id,
                    downloader_job=downloader_job.id,
                    error=str(e),
                )
                downloader_job.success = False
                downloader_job.failure_reason = str(e)
                downloader_job.save()

    def survey(self, source_type=None) -> bool:
        """Retrieves metadata from external source to queue jobs.

        Queries the external source's API to discover an experiment
        and its samples, creates database records for them, and queues
        nomad jobs for them.
        Returns True if successful, False otherwise.
        """
        try:
            experiment, samples = self.discover_experiment_and_samples()
        except Exception:
            logger.exception(
                ("Exception caught while discovering samples. " "Terminating survey job."),
                survey_job=self.survey_job.id,
            )
            self.survey_job.failure_reason = (
                "Exception caught while discovering samples. Terminating survey job."
            )
            return False

        if not experiment:
            logger.info("No experiment found.", survey_job=self.survey_job.id)
            self.survey_job.failure_reason = "No experiment found."
            return False

        try:
            # SRA can have samples with multiple related files,
            # so make sure we download those together.
            if source_type == "SRA":
                for sample in samples:
                    sample_files = sample.original_files.all()
                    self.queue_downloader_job_for_original_files(
                        sample_files, experiment_accession_code=experiment.accession_code
                    )
            else:
                self.queue_downloader_jobs(experiment, samples)
        except Exception:
            logger.exception(
                ("Failed to queue downloader jobs. Terminating survey job."),
                survey_job=self.survey_job.id,
            )
            self.survey_job.failure_reason = (
                "Failed to queue downloader jobs. Terminating survey job."
            )
            return False

        # Update our cached values
        experiment.update_num_samples()
        experiment.update_sample_metadata_fields()
        experiment.update_platform_names()
        experiment.save()

        logger.debug("Survey job completed successfully.", survey_job=self.survey_job.id)
        return True

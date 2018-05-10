import abc
import os

from django.db import transaction
from retrying import retry
from typing import List, Dict

from data_refinery_common.models import (
    DownloaderJob,
    SurveyJob,
    Sample,
    Experiment,
    ExperimentSampleAssociation,
    OriginalFile,
    DownloaderJobOriginalFileAssociation
)
from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class InvalidProcessedFormatError(BaseException):
    pass


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    def downloader_task(self):
        """Returns the Downloaders Enum for the source.
        """
        return Downloaders[self.source_type()]

    @abc.abstractmethod
    def discover_experiments_and_samples(self):
        """Abstract method to survey a source.
        """
        return

    @retry(stop_max_attempt_number=3)
    def queue_downloader_jobs(self, experiment: Experiment):
        """ This enqueues DownloaderJobs on a per-file basis.
        There is a complementary function below foe enqueueing multi-file
        DownloaderJobs.
        """

        # Get all of the undownloaded original files related to this Experiment.
        relations = ExperimentSampleAssociation.objects.filter(experiment=experiment)
        samples = Sample.objects.filter(id__in=relations.values('sample_id'))
        files_to_download = OriginalFile.objects.filter(samples__in=samples.values('pk'), is_downloaded=False)

        downloaded_urls = []
        for original_file in files_to_download:

            # We don't need to create multiple downloaders for the same file.
            if original_file.source_url in downloaded_urls:
              continue

            # There is already a downloader job associated with this file.
            old_assocs = DownloaderJobOriginalFileAssociation.objects.filter(original_file__source_url=original_file.source_url)
            if len(old_assocs) > 0:
              logger.debug("We found an existing DownloaderJob for this file/url.", original_file_id=original_file.id)
              continue

            with transaction.atomic():

              downloader_job = DownloaderJob()
              downloader_job.downloader_task = self.downloader_task()
              downloader_job.accession_code = experiment.accession_code
              downloader_job.save()

              asoc = DownloaderJobOriginalFileAssociation()
              asoc.downloader_job = downloader_job
              asoc.original_file = original_file
              asoc.save()

              downloaded_urls.append(original_file.source_url)

            try:
                logger.info("Queuing downloader job for URL: " + original_file.source_url,
                        survey_job=self.survey_job.id,
                        downloader_job=downloader_job.id)
                send_job(downloader_job.downloader_task, downloader_job.id)
            except:
                # If the task doesn't get sent we don't want the
                # downloader_job to be left floating
                logger.info("Failed to enqueue downloader job for URL: " + original_file.source_url,
                        survey_job=self.survey_job.id,
                        downloader_job=downloader_job.id)
                downloader_job.delete()
                raise

    @retry(stop_max_attempt_number=3)
    def queue_downloader_job_for_original_files(self, original_files: List[OriginalFile]):
        """ Creates a single DownloaderJob with multiple files to download."""

        downloader_job = DownloaderJob()
        downloader_job.downloader_task = self.downloader_task()
        downloader_job.save()

        downloaded_urls = []
        for original_file in original_files:

            # We don't need to create multiple downloaders for the same file.
            if original_file.source_url in downloaded_urls:
              continue

            with transaction.atomic():

              asoc = DownloaderJobOriginalFileAssociation()
              asoc.downloader_job = downloader_job
              asoc.original_file = original_file
              asoc.save()

              downloaded_urls.append(original_file.source_url)

        try:
            logger.info("Queuing downloader job.",
                    survey_job=self.survey_job.id,
                    downloader_job=downloader_job.id)
            send_job(downloader_job.downloader_task, downloader_job.id)
        except:
            # If the task doesn't get sent we don't want the
            # downloader_job to be left floating
            logger.info("Failed to enqueue downloader job.",
                    survey_job=self.survey_job.id,
                    downloader_job=downloader_job.id)
            downloader_job.delete()
            raise

    def survey(self) -> bool:
        try:
            experiment, samples = self.discover_experiment_and_samples()
        except Exception:
            logger.exception(("Exception caught while discovering samples. "
                              "Terminating survey job."),
                             survey_job=self.survey_job.id)
            return False

        if not experiment:
          logger.info("No experiment found.",
                      survey_job=self.survey_job.id)
          return False

        try:
            self.queue_downloader_jobs(experiment)
        except Exception:
            logger.exception(("Failed to queue downloader jobs. "
                              "Terminating survey job."),
                             survey_job=self.survey_job.id)
            return False

        return True

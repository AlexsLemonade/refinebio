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
from data_refinery_common import message_queue, job_lookup, logging


logger = logging.get_and_configure_logger(__name__)


class InvalidProcessedFormatError(Exception):
    pass


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    @abc.abstractmethod
    def discover_experiments_and_samples(self):
        """Abstract method to survey a source.
        """
        return

    @retry(stop_max_attempt_number=3)
    def queue_downloader_jobs(self, experiment: Experiment):
        """ This enqueues DownloaderJobs on a per-file basis.
        There is a complementary function below for enqueueing multi-file
        DownloaderJobs.
        """

        # Get all of the undownloaded original files related to this Experiment.
        relations = ExperimentSampleAssociation.objects.filter(experiment=experiment)
        samples = Sample.objects.filter(id__in=relations.values('sample_id'))
        files_to_download = OriginalFile.objects.filter(
            samples__in=samples.values('pk'), is_downloaded=False)

        downloaded_urls = []
        for original_file in files_to_download:

            # We don't need to create multiple downloaders for the same file.
            if original_file.source_url in downloaded_urls:
                continue

            # There is already a downloader job associated with this file.
            old_assocs = DownloaderJobOriginalFileAssociation.objects.filter(
                original_file__source_url=original_file.source_url)
            if len(old_assocs) > 0:
                logger.debug("We found an existing DownloaderJob for this file/url.",
                             original_file_id=original_file.id)
                continue

            sample_object = original_file.samples.first()
            downloader_task = job_lookup.determine_downloader_task(sample_object)

            if downloader_task == job_lookup.Downloaders.NONE:
                logger.info("No valid downloader task found for sample.",
                            sample=sample_object.id,
                            original_file=original_file.id)
            else:
                with transaction.atomic():
                    downloader_job = DownloaderJob()
                    downloader_job.downloader_task = downloader_task.value
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
                    message_queue.send_job(downloader_task, downloader_job.id)
                except Exception as e:
                    # If the task doesn't get sent we don't want the
                    # downloader_job to be left floating
                    logger.exception("Failed to enqueue downloader job for URL: "
                                     + original_file.source_url,
                                     survey_job=self.survey_job.id,
                                     downloader_job=downloader_job.id)
                    downloader_job.success = False
                    downloader_job.failure_reason = str(e)
                    downloader_job.save()

    @retry(stop_max_attempt_number=3)
    def queue_downloader_job_for_original_files(self,
                                                original_files: List[OriginalFile],
                                                experiment_accession_code: str=None,
                                                is_transcriptome: bool=False
                                                ):
        """ Creates a single DownloaderJob with multiple files to download.
        """
        # Transcriptome is a special case because there's no sample_object.
        if is_transcriptome:
            downloader_task = job_lookup.Downloaders.TRANSCRIPTOME_INDEX
        else:
            sample_object = original_files[0].samples.first()
            downloader_task = job_lookup.determine_downloader_task(sample_object)

        if downloader_task == job_lookup.Downloaders.NONE:
            logger.info("No valid downloader task found for sample.",
                        sample=sample_object.id,
                        original_file=original_files[0].id)
        else:
            downloader_job = DownloaderJob()
            downloader_job.downloader_task = downloader_task.value
            downloader_job.accession_code = experiment_accession_code
            downloader_job.save()

            downloaded_urls = []
            for original_file in original_files:
                asoc = DownloaderJobOriginalFileAssociation()
                asoc.downloader_job = downloader_job
                asoc.original_file = original_file
                asoc.save()

                downloaded_urls.append(original_file.source_url)

            try:
                logger.info("Queuing downloader job.",
                            survey_job=self.survey_job.id,
                            downloader_job=downloader_job.id,
                            downloaded_urls=downloaded_urls)
                message_queue.send_job(downloader_task, downloader_job)
            except Exception as e:
                # If the task doesn't get sent we don't want the
                # downloader_job to be left floating
                logger.exception("Failed to enqueue downloader job.",
                                 survey_job=self.survey_job.id,
                                 downloader_job=downloader_job.id)
                downloader_job.success = False
                downloader_job.failure_reason = str(e)
                downloader_job.save()

    def survey(self, source_type=None) -> bool:
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

            # SRA can have samples with multiple related files,
            # so make sure we download those together.
            if source_type == "SRA":
                for sample in samples:
                    sample_files = sample.original_files.all()
                    self.queue_downloader_job_for_original_files(sample_files,
                                                                 experiment_accession_code=experiment.accession_code)
            else:
                self.queue_downloader_jobs(experiment)
        except Exception:
            logger.exception(("Failed to queue downloader jobs. "
                              "Terminating survey job."),
                             survey_job=self.survey_job.id)
            return False

        logger.info("Survey job completed successfully." survey_job=self.survey_job.id)
        return True

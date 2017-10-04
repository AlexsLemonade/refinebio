import abc
from typing import List, Dict
from retrying import retry
from django.db import transaction
from data_refinery_common.models import (
    Batch,
    BatchKeyValue,
    BatchStatuses,
    File,
    DownloaderJob,
    SurveyJob
)
from data_refinery_foreman.surveyor.message_queue import app
from data_refinery_common.job_lookup import DOWNLOADER_TASK_LOOKUP


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class InvalidProcessedFormatError(BaseException):
    pass


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta
    batches = []

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    def group_batches(self) -> List[List[Batch]]:
        """Groups batches together which should be downloaded together.

        The default implementation groups batches based on the
        download URL of their first File.
        """
        groups = []
        download_urls = {batch.files[0].download_url for batch in self.batches}
        for url in download_urls:
            groups.append([batch for batch in self.batches if batch.files[0].download_url == url])

        return groups

    @abc.abstractmethod
    def determine_pipeline(self,
                           batch: Batch,
                           files: List[File],
                           key_values: Dict = {}):
        """Determines the appropriate pipeline for the batch.

        Returns a string that represents a processor pipeline.
        Must return a member of PipelineEnums.
        """
        return

    def downloader_task(self):
        """Returns the downloader task for the source.

        Returns the Celery Downloader Task name from the
        data_refinery_workers project which should be queued to
        download Batches discovered by this surveyor.
        """
        return DOWNLOADER_TASK_LOOKUP[self.source_type()]

    @retry(stop_max_attempt_number=3)
    @transaction.atomic
    def add_batch(self,
                  platform_accession_code: str,
                  experiment_accession_code: str,
                  organism_id: int,
                  organism_name: str,
                  experiment_title: str,
                  release_date,
                  last_uploaded_date,
                  files: List[File],
                  size_in_bytes: int = -1,
                  key_values: Dict = {}):
        # Prevent creating duplicate Batches.
        for file in files:
            if File.objects.filter(name=file.name).count() != 0:
                logger.info(("Skipping sample with name %s because a File already exists with"
                             "that name."),
                            file.name)
                return

        batch = Batch(survey_job=self.survey_job,
                      source_type=self.source_type(),
                      status=BatchStatuses.NEW.value,
                      size_in_bytes=-1,
                      platform_accession_code=platform_accession_code,
                      experiment_accession_code=experiment_accession_code,
                      organism_id=organism_id,
                      organism_name=organism_name,
                      experiment_title=experiment_title,
                      release_date=release_date,
                      last_uploaded_date=last_uploaded_date)

        batch.pipeline_required = self.determine_pipeline(batch, files, key_values)
        batch.save()

        for file in files:
            file.batch = batch
            file.save()

        for key, value in key_values.items():
            BatchKeyValue(batch=batch,
                          key=key,
                          value=value).save()

        self.batches.append(batch)

    @abc.abstractmethod
    def discover_batches(self):
        """Abstract method to survey a source.

        Implementations of this method should do the following:
        1. Query the external source to discover batches that should be
           downloaded.
        2. Create a Batch object for each discovered batch and optionally
           a list of BatchKeyValues.
        3. Call self.handle_batch for each Batch object that is created.

        Each Batch object should have the following fields populated:
            size_in_bytes
            download_url
            raw_format -- it is possible this will not yet be known
            accession_code
            organism

        The following fields will be set by handle_batch:
            survey_job
            source_type
            pipeline_required
            status

        The processed_format should be set if it is already known what it
        will be. If it is not set then the determine_pipeline method must
        return a DiscoveryPipeline (a pipeline that determines what the
        raw_format of the data is, what the processed_format should be, and
        which pipeline to use to transform it).

        Return:
        This method should return True if the job completed successfully with
        no issues. Otherwise it should log any issues and return False.
        """
        return

    @retry(stop_max_attempt_number=3)
    def queue_downloader_jobs(self, batches: List[Batch]):
        if len(batches) > 0:
            downloader_task = self.downloader_task()

            with transaction.atomic():
                downloader_job = DownloaderJob.create_job_and_relationships(
                    batches=batches, downloader_task=downloader_task)

            logger.info("Survey job #%d is queuing downloader job #%d.",
                        self.survey_job.id,
                        downloader_job.id)
            try:
                app.send_task(downloader_task, args=[downloader_job.id])
            except:
                # If the task doesn't get sent we don't want the
                # downloader_job to be left floating
                downloader_job.delete()
                raise
        else:
            logger.info("Survey job #%d found no new Batches.",
                        self.survey_job.id)

    def survey(self) -> bool:
        try:
            self.discover_batches()
        except Exception:
            logger.exception(("Exception caught while discovering batches. "
                              "Terminating survey job #%d."),
                             self.survey_job.id)
            return False

        for group in self.group_batches():
            try:
                self.queue_downloader_jobs()
            except Exception:
                logger.exception(("Failed to queue downloader jobs. "
                                  "Terminating survey job #%d."),
                                 self.survey_job.id)
                return False

        return True

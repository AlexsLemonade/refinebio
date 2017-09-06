import abc
import os
from typing import List
from retrying import retry
from django.db import transaction
from data_refinery_models.models import (
    Batch,
    BatchStatuses,
    DownloaderJob,
    SurveyJob
)
from data_refinery_foreman.surveyor.message_queue import app
from data_refinery_common.job_lookup import DiscoveryPipeline, DOWNLOADER_TASK_LOOKUP


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class InvalidProcessedFormatError(BaseException):
    pass


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    @abc.abstractmethod
    def determine_pipeline(self,
                           batch: Batch):
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

    def handle_batches(self, batches: List[Batch]):
        new_batches = []
        for batch in batches:
            if not Batch.is_new_batch(batch):
                logger.info(("Skipping sample with name %s because a batch already exists with"
                             "that name."),
                            batch.name)
                continue

            batch.survey_job = self.survey_job
            batch.source_type = self.source_type()
            batch.status = BatchStatuses.NEW.value

            pipeline_required = self.determine_pipeline(batch)
            if (pipeline_required is DiscoveryPipeline) or batch.processed_format:
                batch.pipeline_required = pipeline_required.value
            else:
                message = ("Batches must have the processed_format field set "
                           "unless the pipeline returned by determine_pipeline "
                           "is of the type DiscoveryPipeline.")
                raise InvalidProcessedFormatError(message)

            batch.internal_location = os.path.join(batch.platform_accession_code,
                                                   batch.pipeline_required)

            batch.save()
            new_batches.append(batch)

        @retry(stop_max_attempt_number=3)
        def save_batches_start_job():
            if len(new_batches) > 0:
                downloader_task = self.downloader_task()

                with transaction.atomic():
                    downloader_job = DownloaderJob.create_job_and_relationships(
                        batches=new_batches, downloader_task=downloader_task)

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
                logger.info("Survey job #% found no new Batches.",
                            self.survey_job.id)

        try:
            save_batches_start_job()
        except Exception:
            logger.exception(("Failed to save batches to database three times. "
                              "Terminating survey job #%d."),
                             self.survey_job.id)
            raise

    @abc.abstractmethod
    def survey(self):
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

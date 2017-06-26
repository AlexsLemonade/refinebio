import abc
import os
from enum import Enum
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


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class InvalidProcessedFormatError(BaseException):
    pass


class PipelineEnums(Enum):
    """An abstract class to enumerate valid processor pipelines.

    Enumerations which extend this class are valid values for the
    pipeline_required field of the Batches table.
    """
    pass


class ProcessorPipeline(PipelineEnums):
    """Pipelines which perform some kind of processing on the data."""
    AFFY_TO_PCL = "AFFY_TO_PCL"


class DiscoveryPipeline(PipelineEnums):
    """Pipelines which discover appropriate processing for the data."""
    pass


class ExternalSourceSurveyor:
    __metaclass__ = abc.ABCMeta

    def __init__(self, survey_job: SurveyJob):
        self.survey_job = survey_job

    @abc.abstractproperty
    def source_type(self):
        return

    @abc.abstractproperty
    def downloader_task(self):
        """Abstract property representing the downloader task.

        Should return the Celery Downloader Task name from the
        data_refinery_workers project which should be queued to
        download Batches discovered by this surveyor.
        """
        return

    @abc.abstractmethod
    def determine_pipeline(self,
                           batch: Batch):
        """Determines the appropriate pipeline for the batch.

        Returns a string that represents a processor pipeline.
        Must return a member of PipelineEnums.
        """
        return

    def handle_batches(self, batches: List[Batch]):
        for batch in batches:
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

        @retry(stop_max_attempt_number=3)
        @transaction.atomic
        def save_batches_start_job():
            downloader_task = self.downloader_task()
            downloader_job = DownloaderJob.create_job_and_relationships(
                batches=batches, downloader_task=downloader_task)
            app.send_task(downloader_task, args=[downloader_job.id])

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

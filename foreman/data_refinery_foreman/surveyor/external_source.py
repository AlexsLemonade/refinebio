import abc
from enum import Enum
from typing import List
from data_refinery_models.models import (
    Batch,
    BatchStatuses,
    BatchKeyValue,
    SurveyJob
)


class PipelineEnums(Enum):
    """An abstract class to enumerate valid processor pipelines.

    Enumerations which extend this class are valid values for the
    pipeline_required field of the Batches table."""
    pass


class ProcessorPipeline(PipelineEnums):
    """Pipelines which perform some kind of processing on the data."""
    MICRO_ARRAY_TO_PCL = "MICRO_ARRAY_TO_PCL"


class DiscoveryPipeline(PipelineEnums):
    """Pipelines which discover what kind of processing is appropriate
    for the data."""
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
                           batch: Batch,
                           key_values: List[BatchKeyValue] = []):
        """Determines the appropriate processor pipeline for the batch
        and returns a string that represents a processor pipeline.
        Must return a member of PipelineEnums."""
        return

    def handle_batch(self, batch: Batch, key_values: BatchKeyValue = None):
        batch.survey_job = self.survey_job
        batch.source_type = self.source_type()
        batch.status = BatchStatuses.NEW.value
        batch.internal_location = (batch.accession_code + "/"
                                   + batch.pipeline_required + "/")

        pipeline_required = self.determine_pipeline(batch, key_values)
        if (pipeline_required is DiscoveryPipeline) or batch.processed_format:
            batch.pipeline_required = pipeline_required.value
        else:
            message = ("Batches must have the processed_format field set " +
                       "unless the pipeline returned by determine_pipeline" +
                       "is of the type DiscoveryPipeline.")
            # Also should be more specific
            raise Exception(message)

        # This is also where we will queue the downloader job
        if batch.save():
            return True
        else:
            return False

    @abc.abstractmethod
    def survey(self, survey_job: SurveyJob):
        """Implementations of this function should do the following:
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

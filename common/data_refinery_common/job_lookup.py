from enum import Enum


class PipelineEnums(Enum):
    """An abstract class to enumerate valid processor pipelines.

    Enumerations which extend this class are valid values for the
    pipeline_required field of the Batches table.
    """
    pass


class ProcessorPipeline(PipelineEnums):
    """Pipelines which perform some kind of processing on the data."""
    AFFY_TO_PCL = "AFFY_TO_PCL"
    AGILENT_TWOCOLOR_TO_PCL = "AGILENT_TWOCOLOR_TO_PCL"
    SALMON = "SALMON"
    ILLUMINA_TO_PCL = "ILLUMINA_TO_PCL"
    TRANSCRIPTOME_INDEX_LONG = "TRANSCRIPTOME_INDEX_LONG"
    TRANSCRIPTOME_INDEX_SHORT = "TRANSCRIPTOME_INDEX_SHORT"
    NO_OP = "NO_OP"
    NONE = "NONE"


class DiscoveryPipeline(PipelineEnums):
    """Pipelines which discover appropriate processing for the data."""
    pass


class Downloaders(Enum):
    """An enumeration of downloaders for downloader_task."""
    ARRAY_EXPRESS = "ARRAY_EXPRESS"
    SRA = "SRA"
    TRANSCRIPTOME_INDEX = "TRANSCRIPTOME_INDEX"
    GEO = "GEO"

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
    SALMON = "SALMON"
    TRANSCRIPTOME_INDEX = "TRANSCRIPTOME_INDEX"
    NO_OP = "NO_OP"
    NONE = "NONE"


class DiscoveryPipeline(PipelineEnums):
    """Pipelines which discover appropriate processing for the data."""
    pass


# Maps processor pipeline names to the Celery task definitions.
PROCESSOR_PIPELINE_LOOKUP = {
    "AFFY_TO_PCL": "data_refinery_workers.processors.array_express.affy_to_pcl",
    "NO_OP": "data_refinery_workers.processors.no_op.no_op_processor",
    "SALMON": "data_refinery_workers.processors.salmon.salmon",
    "TRANSCRIPTOME_INDEX": "data_refinery_workers.processors.transcriptome_index.build_index",
}


class Downloaders(Enum):
    """An enumeration of downloaders for batch.downloader_task."""
    ARRAY_EXPRESS = "ARRAY_EXPRESS"
    SRA = "SRA"
    TRANSCRIPTOME_INDEX = "TRANSCRIPTOME_INDEX"


DOWNLOADER_TASK_LOOKUP = {
    "ARRAY_EXPRESS": "data_refinery_workers.downloaders.array_express.download_array_express",
    "SRA": "data_refinery_workers.downloaders.sra.download_sra",
    "TRANSCRIPTOME_INDEX": ("data_refinery_workers.downloaders."
                            "transcriptome_index.download_transcriptome"),
}

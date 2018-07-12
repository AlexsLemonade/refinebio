"""This module defines the names of pipelines and processors."""

from enum import Enum, unique

@unique
class PipelineEnum(Enum):
    AGILENT_TWOCOLOR = "Agilent Two Color"
    ARRAY_EXPRESS = "Array Express"
    ILLUMINA = "Illumina"
    NO_OP = "No Op"
    SALMON = "Salmon"
    SMASHER = "Smasher"
    TX_INDEX = "Transcriptome Index"


@unique
class ProcessorEnum(Enum):
    # One processor in "Agilent Two Color" pipeline:
    AGILENT_TWOCOLOR = "Agilent SCAN TwoColor"

    # One processor in "Array Express" pipeline:
    AFFYMETRIX_SCAN = "Affymetrix SCAN"

    # One processor in "Illumina" pipeline:
    ILLUMINA_SCAN = "Illumina SCAN"

    # One processor in "No Op" pipeline:
    SUBMITTER_PROCESSED = "Submitter-processed"

    # Four processors in "Salmon" pipeline:
    TXIMPORT = "Tximport"
    SALMON_QUANT = "Salmon Quant"
    MULTIQC = "MultiQC"
    SALMONTOOLS = "Salmontools"

    # No processors in "Smasher" pipeline (yet)

    # One processor in "Transcriptome Index" pipeline:
    TX_INDEX = "Transcriptome Index"

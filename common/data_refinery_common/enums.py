from enum import Enum, unique


@unique
class PipelineEnum(Enum):
    """Hardcoded pipeline names."""

    AGILENT_TWOCOLOR = "Agilent Two Color"
    ARRAY_EXPRESS = "Array Express"
    ILLUMINA = "Illumina"
    NO_OP = "No Op"
    SALMON = "Salmon"
    SMASHER = "Smasher"
    TXIMPORT = "Tximport"
    TX_INDEX = "Transcriptome Index"
    QN_REFERENCE = "Quantile Normalization Reference"
    JANITOR = "Janitor"
    CREATE_COMPENDIA = "Compendia"
    CREATE_QUANTPENDIA = "Quantpendia"


@unique
class ProcessorEnum(Enum):
    """Hardcoded processor info in each pipeline."""

    # One processor in "Agilent Two Color" pipeline
    AGILENT_TWOCOLOR = {
        "name": "Agilent SCAN TwoColor",
        "docker_img": "not_available_yet",
        "yml_file": "agilent_twocolor.yml",
    }

    # One processor in "Array Express" pipeline
    AFFYMETRIX_SCAN = {
        "name": "Affymetrix SCAN",
        "docker_img": "dr_affymetrix",
        "yml_file": "affymetrix.yml",
    }

    # One processor in "Illumina" pipeline
    ILLUMINA_SCAN = {
        "name": "Illumina SCAN",
        "docker_img": "dr_illumina",
        "yml_file": "illumina.yml",
    }

    # One processor in "No Op" pipeline
    SUBMITTER_PROCESSED = {
        "name": "Submitter-processed",
        "docker_img": "dr_no_op",
        "yml_file": "no_op.yml",
    }

    # Three processors in "Salmon" pipeline
    SALMON_QUANT = {
        "name": "Salmon Quant",
        "docker_img": "dr_salmon",
        "yml_file": "salmon_quant.yml",
    }
    SALMONTOOLS = {"name": "Salmontools", "docker_img": "dr_salmon", "yml_file": "salmontools.yml"}
    TXIMPORT = {"name": "Tximport", "docker_img": "dr_salmon", "yml_file": "tximport.yml"}

    # One processor in "Smasher" pipeline
    SMASHER = {"name": "Smasher", "docker_img": "dr_smasher", "yml_file": "smasher.yml"}

    # One processor in "Transcriptome Index" pipeline
    TX_INDEX = {
        "name": "Transcriptome Index",
        "docker_img": "dr_transcriptome",
        "yml_file": "transcriptome_index.yml",
    }

    QN_REFERENCE = {
        "name": "Quantile Normalization Reference",
        "docker_img": "dr_smasher",
        "yml_file": "qn.yml",
    }

    CREATE_COMPENDIA = {
        "name": "Compendia Creation",
        "docker_img": "dr_compendia",
        "yml_file": "compendia.yml",
    }

    CREATE_QUANTPENDIA = {
        "name": "Quantpendia Creation",
        "docker_img": "dr_compendia",
        "yml_file": "compendia.yml",
    }

    @classmethod
    def has_key(cls, key):
        """Class method that tells whether a certain key exists."""
        return key in cls.__members__


class ProcessorPipeline(Enum):
    """An enumeration of supported processors"""

    AFFY_TO_PCL = "AFFY_TO_PCL"
    AGILENT_ONECOLOR_TO_PCL = "AGILENT_ONECOLOR_TO_PCL"  # Currently unsupported
    AGILENT_TWOCOLOR_TO_PCL = "AGILENT_TWOCOLOR_TO_PCL"
    SALMON = "SALMON"
    TXIMPORT = "TXIMPORT"
    ILLUMINA_TO_PCL = "ILLUMINA_TO_PCL"
    TRANSCRIPTOME_INDEX_LONG = "TRANSCRIPTOME_INDEX_LONG"
    TRANSCRIPTOME_INDEX_SHORT = "TRANSCRIPTOME_INDEX_SHORT"
    SMASHER = "SMASHER"
    NO_OP = "NO_OP"
    QN_REFERENCE = "QN_REFERENCE"
    JANITOR = "JANITOR"
    NONE = "NONE"
    CREATE_COMPENDIA = "CREATE_COMPENDIA"
    CREATE_QUANTPENDIA = "CREATE_QUANTPENDIA"


SMASHER_JOB_TYPES = [
    ProcessorPipeline.SMASHER,
    ProcessorPipeline.QN_REFERENCE,
    ProcessorPipeline.CREATE_COMPENDIA,
    ProcessorPipeline.CREATE_QUANTPENDIA,
]


class Downloaders(Enum):
    """An enumeration of downloaders for downloader_task."""

    ARRAY_EXPRESS = "ARRAY_EXPRESS"
    SRA = "SRA"
    TRANSCRIPTOME_INDEX = "TRANSCRIPTOME_INDEX"
    GEO = "GEO"
    NONE = "NONE"


class SurveyJobTypes(Enum):
    """An enumeration of downloaders for downloader_task."""

    SURVEYOR = "SURVEYOR"

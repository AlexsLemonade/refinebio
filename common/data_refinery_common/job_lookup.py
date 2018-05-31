from enum import Enum
from typing import List

from data_refinery_common import utils
from data_refinery_common.models import OriginalFile


class PipelineEnums(Enum):

    """An abstract class to enumerate valid processor pipelines.

    Enumerations which extend this class are valid values for the
    pipeline_required field of the Batches table.
    """
    pass


class ProcessorPipeline(PipelineEnums):

    """An enumeration of supported processors"""
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
    NONE = "NONE"


def _is_platform_supported(platform: str) -> bool:
    upper_platform = platform.upper()

    # Check if this is a supported Microarray platform.
    for supported_platform in utils.get_supported_microarray_platforms():
        if (supported_platform["platform_accession"].upper() == upper_platform
                or supported_platform["external_accession"].upper() == upper_platform):
            return True

    # Check if this is a supported RNASeq platform.
    # GEO RNASeq platform titles often have organisms appended to
    # an otherwise recognizable platform. The list of supported
    # RNASeq platforms isn't long, so see if any of them are
    # contained within what GEO gave us.
    # Example: GSE69572 has a platform title of:
    # 'Illumina Genome Analyzer IIx (Glycine max)'
    # Which should match 'Illumina Genome Analyzer IIx'
    # because RNASeq platforms are organism agnostic.
    for supported_platform in utils.get_supported_rnaseq_platforms():
        if supported_platform.upper() in upper_platform:
            return True

    return False


def _determine_microarray_manufacturer(platform: str) -> str:
    upper_platform = platform.upper()
    for supported_platform in utils.get_supported_microarray_platforms():
        # Check both the internal and external accessions so we don't
        # have to worry about the wrong accession getting passed in.
        upper_internal_accession = supported_platform["platform_accession"].upper()
        if (platform == upper_internal_accession
                or platform == supported_platform["external_accession"].upper()):
            # Found the right accession, now determine if the manufacturer is Illumina:
            if "ILLUMINA" in supported_platform["platform_accession"]:
                return "ILLUMINA"
            else:
                # It's not Illumina, the only other supported Microarray platform is Affy:
                return "AFFYMETRIX"

    # This is easily prevented by only calling this function after _is_platform_supported()
    raise UnsupportedPlatformError(
        "Cannot determine manufacterer for unsupported platform: " + platform)


def determine_downloader_task(source_database: str,
                              technology: str,
                              has_raw: bool,
                              platform: str,
                              original_files: List[OriginalFile]=[]) -> Downloaders:
    if _is_platform_supported(platform):
        return Downloaders[source_database]
    elif has_raw:
        for original_file in original_files:
            if original_file.source_filename[-4:].upper() == ".CEL":
                return Downloaders[source_database]

    return Downloaders.NONE


def determine_processor_pipeline(source_database: str,
                                 technology: str,
                                 has_raw: bool,
                                 platform: str) -> ProcessorPipeline:
    if not _is_platform_supported(platform):
        return ProcessorPipeline.NONE

    if technology == "MICROARRAY":
        if not has_raw:
            return ProcessorPipeline.NO_OP
        else:
            try:
                manufacterer = _determine_microarray_manufacturer(platform)
            except UnsupportedPlatformError:
                # This shouldn't happen since we first check that the
                # platform is supported.
                return ProcessorPipeline.NONE

            if manufacterer == "ILLUMINA":
                return ProcessorPipeline.ILLUMINA_TO_PCL
            elif manufacterer == "AFFYMETRIX":
                return ProcessorPipeline.AFFY_TO_PCL

    else:
        if not has_raw:
            # Only NO_OP RNASeq data if it's from SRA, because
            # ArrayExpress and GEO often only have processed data
            # while linking to SRA for the raw data.
            if source_database == "SRA":
                return ProcessorPipeline.NO_OP
            else:
                return ProcessorPipeline.NONE
        else:
            return ProcessorPipeline.SALMON

    # Shouldn't get here, but just in case
    return ProcessorPipeline.NONE

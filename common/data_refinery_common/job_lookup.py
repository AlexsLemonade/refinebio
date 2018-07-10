from enum import Enum
from typing import List

from data_refinery_common import utils
from data_refinery_common.models import Sample, OriginalFile, OriginalFileSampleAssociation
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


class PipelineEnums(Enum):

    """An abstract class to enumerate valid processor pipelines.

    Enumerations which extend this class are valid values for the
    pipeline_required field of the Batches table.
    """
    pass


class ProcessorPipeline(PipelineEnums):

    """An enumeration of supported processors"""
    AFFY_TO_PCL = "AFFY_TO_PCL"
    AGILENT_ONECOLOR_TO_PCL = "AGILENT_ONECOLOR_TO_PCL"  # Currently unsupported
    AGILENT_TWOCOLOR_TO_PCL = "AGILENT_TWOCOLOR_TO_PCL"
    SALMON = "SALMON"
    ILLUMINA_TO_PCL = "ILLUMINA_TO_PCL"
    TRANSCRIPTOME_INDEX_LONG = "TRANSCRIPTOME_INDEX_LONG"
    TRANSCRIPTOME_INDEX_SHORT = "TRANSCRIPTOME_INDEX_SHORT"
    SMASHER = "SMASHER"
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
    """Determines if platform is a platform_accession we support or not.

    It does so by trying to correct for common string issues such as
    case and spacing and then comparing against our configuration
    files which specify which platform are supported.
    """
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
        # Spacing can be inconsistent, easiest to just remove it entirely.
        if supported_platform.upper().replace(" ", "") in upper_platform.replace(" ", ""):
            return True

    return False


def determine_downloader_task(sample_object: Sample) -> Downloaders:
    """Returns the Downloaders enum appropriate for the sample.

    For any sample which has a supported platform this is entirely
    based on the source database since we have a one-to-one mapping
    between sources and downloader tasks. If the platform isn't
    supported, then we don't want to download the sample so we return
    Downloaders.NONE. However any sample with a .CEl file could have
    inaccurate platform information which we potentially correct for
    after downloading it.
    """
    if _is_platform_supported(sample_object.platform_accession_code):
        return Downloaders[sample_object.source_database]
    elif sample_object.has_raw:
        # Sometimes Array Express lies about what a sample's platform
        # is. Therefore, if there's a .CEL file we'll download it and
        # determine the platform from that.
        relations = OriginalFileSampleAssociation.objects.filter(sample=sample_object)
        original_files = OriginalFile.objects.filter(id__in=relations.values('original_file_id'))
        for original_file in original_files:
            if original_file.source_filename[-4:].upper() == ".CEL":
                return Downloaders[sample_object.source_database]

    return Downloaders.NONE


def determine_processor_pipeline(sample_object: Sample) -> ProcessorPipeline:
    """Determines the appropriate processor pipeline for the sample.

    This is mostly a giant set of nested if statements, so describing
    the logic wouldn't add very much. However, the general flow is:
      - Is the platform supported? If not return NONE cause we don't want it.
      - Does it have raw data? If not NO_OP the data.
        (With one exception explained in comments.)
      - Is it Microarray data? If so determine which processor based on its
        manufacturer Otherwise it's SALMON-time.
    """
    if not _is_platform_supported(sample_object.platform_accession_code):
        return ProcessorPipeline.NONE

    if sample_object.technology == "MICROARRAY":
        if not sample_object.has_raw:
            return ProcessorPipeline.NO_OP
        else:
            if sample_object.manufacturer == "ILLUMINA":
                return ProcessorPipeline.ILLUMINA_TO_PCL
            elif sample_object.manufacturer == "AFFYMETRIX":
                return ProcessorPipeline.AFFY_TO_PCL
            elif sample_object.manufacturer == "AGILENT":
                # We currently aren't prepared to process Agilent because we don't have
                # whitelist of supported platforms for it. However this code works so
                # let's keep it around until we're ready for Agilent.
                annotations = sample_object.sampleannotation_set.all()[0]
                channel1_protocol = annotations.data.get('label_protocol_ch1', "").upper()
                channel2_protocol = annotations.data.get('label_protocol_ch2', "").upper()
                if ('AGILENT' in channel1_protocol) and ('AGILENT' in channel2_protocol):
                    return ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL
                else:
                    return ProcessorPipeline.AGILENT_ONECOLOR_TO_PCL
            elif sample_object.manufacturer == "UNKNOWN":
                logger.error("Found a Sample on a supported platform with an unknown manufacturer.",
                             sample=sample_object.id,
                             platform_accession=sample_object.platform_accession_code,
                             accession=sample_object.accession_code,
                             manufacturer=sample_object.manufacturer,
                             platform_name=sample_object.platform_name
                            )
                return ProcessorPipeline.NONE

    elif sample_object.has_raw:
        return ProcessorPipeline.SALMON

    # Shouldn't get here, but just in case
    return ProcessorPipeline.NONE

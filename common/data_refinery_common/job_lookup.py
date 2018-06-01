from enum import Enum
from typing import List

from data_refinery_common import utils
from data_refinery_common.models import OriginalFile
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


def determine_downloader_task(sample_object: Sample) -> Downloaders:
    if _is_platform_supported(sample_object.platform_accession_code):
        return Downloaders[source_database]
    elif sample_object.has_raw:
        # Sometimes Array Express lies about what a sample's platform
        # is. Therefore, if there's a .CEL file we'll download it and
        # determine the platform from that.
        relations = OriginalFileSampleAssociation.objects.filter(sample=sample_object)
        original_files = OriginalFile.objects.filter(id__in=relations.values('original_file_id'))
        for original_file in original_files:
            if original_file.source_filename[-4:].upper() == ".CEL":
                return Downloaders[source_database]

    return Downloaders.NONE


def determine_processor_pipeline(sample_object: Sample) -> ProcessorPipeline:
    if not _is_platform_supported(sample_object.platform_accession_code):
        return ProcessorPipeline.NONE

    if sample_object.technology == "MICROARRAY":
        if not sample_object.has_raw:
            return ProcessorPipeline.NO_OP
        else:
            try:
                manufacterer = _determine_microarray_manufacturer(
                    sample_object.platform_accession_code)
            except UnsupportedPlatformError:
                # This shouldn't happen since we first check that the
                # platform is supported.
                return ProcessorPipeline.NONE

            if sample_object.manufacterer == "ILLUMINA":
                return ProcessorPipeline.ILLUMINA_TO_PCL
            elif sample_object.manufacterer == "AFFYMETRIX":
                return ProcessorPipeline.AFFY_TO_PCL
            elif sample_object.manufacterer == "AGILENT":
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
            elif sample_object.manufacterer == "UNKNOWN":
                logger.error("Found a Sample on a supported platform with an unknown manufacterer."
                             sample=sample_object.id,
                             platform_accession=sample_object.platform_accession_code)
                return ProcessorPipeline.NONE

    else:
        if not sample_object.has_raw:
            # Only NO_OP RNASeq data if it's from SRA, because
            # ArrayExpress and GEO often only have processed data
            # while linking to SRA for the raw data.
            if sample_object.source_database == "SRA":
                return ProcessorPipeline.NO_OP
            else:
                return ProcessorPipeline.NONE
        else:
            return ProcessorPipeline.SALMON

    # Shouldn't get here, but just in case
    return ProcessorPipeline.NONE

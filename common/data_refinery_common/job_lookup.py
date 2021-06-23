import re

from data_refinery_common import utils
from data_refinery_common.enums import Downloaders, ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ProcessorJob, Sample

logger = get_and_configure_logger(__name__)


def does_processor_job_have_samples(job: ProcessorJob):
    return not (
        job.pipeline_applied == ProcessorPipeline.SMASHER.value
        or job.pipeline_applied == ProcessorPipeline.JANITOR.value
        or job.pipeline_applied == ProcessorPipeline.QN_REFERENCE.value
    )


def is_file_rnaseq(filename: str) -> bool:
    """Returns true if `filename` matches the pattern of an RNAseq file, false otherwise."""
    if not filename:
        return False

    matches_accession_pattern = re.search("[SED]RR\d{5,}", filename)

    return (
        filename[-5:].upper() == "FASTQ"
        or filename[-8:].upper() == "FASTQ.GZ"
        or filename[-2:].upper() == "FQ"
        or filename[-3:].upper() == "SRA"
        or filename[-5:].upper() == "FQ.GZ"
        or matches_accession_pattern
    )


def _is_platform_supported(platform: str) -> bool:
    """Determines if platform is a platform_accession we support or not.

    It does so by trying to correct for common string issues such as
    case and spacing and then comparing against our configuration
    files which specify which platform are supported.
    """
    upper_platform = platform.upper()

    # Check if this is a supported Microarray platform.
    for supported_platform in utils.get_supported_microarray_platforms():
        if (
            supported_platform["platform_accession"].upper() == upper_platform
            or supported_platform["external_accession"].upper() == upper_platform
        ):
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
        has_cel_original_files = sample_object.original_files.filter(
            source_filename__iendswith=".CEL"
        ).exists()
        if has_cel_original_files:
            return Downloaders[sample_object.source_database]

    return Downloaders.NONE


def determine_processor_pipeline(sample_object: Sample, original_file=None) -> ProcessorPipeline:
    """Determines the appropriate processor pipeline for the sample.

    This is mostly a giant set of nested if statements, so describing
    the logic wouldn't add very much. However, the general flow is:
      - Is its file extension .CEL or FQ*/FASTQ*? Send it to the correct processor.
      - Does it have raw data? If not NO_OP the data.
        (With one exception explained in comments.)
      - Is the platform supported? If not return NONE cause we don't want it.
      - Is it Microarray data? If so determine which processor based on its
        manufacturer Otherwise it's SALMON-time.
    """
    if original_file:
        if original_file.is_affy_data():
            return ProcessorPipeline.AFFY_TO_PCL
        if is_file_rnaseq(original_file.filename):
            return ProcessorPipeline.SALMON

    # We NO_OP processed data. It's what we do.
    if not sample_object.has_raw or (original_file and ".processed" in original_file.source_url):
        return ProcessorPipeline.NO_OP

    if not _is_platform_supported(sample_object.platform_accession_code):
        return ProcessorPipeline.NONE

    if sample_object.technology == "MICROARRAY":
        if sample_object.manufacturer == "ILLUMINA":
            return ProcessorPipeline.ILLUMINA_TO_PCL
        elif sample_object.manufacturer == "AFFYMETRIX":
            # Optional explicit filetype checks
            if original_file and original_file.filename[-3:].upper() == "TXT":
                return ProcessorPipeline.NO_OP
            else:
                return ProcessorPipeline.AFFY_TO_PCL
        elif sample_object.manufacturer == "AGILENT":
            # We currently aren't prepared to process Agilent because we don't have
            # whitelist of supported platforms for it. However this code works so
            # let's keep it around until we're ready for Agilent.
            annotations = sample_object.sampleannotation_set.all()[0]
            channel1_protocol = annotations.data.get("label_protocol_ch1", "").upper()
            channel2_protocol = annotations.data.get("label_protocol_ch2", "").upper()
            if ("AGILENT" in channel1_protocol) and ("AGILENT" in channel2_protocol):
                return ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL
            else:
                return ProcessorPipeline.AGILENT_ONECOLOR_TO_PCL
        elif sample_object.manufacturer == "UNKNOWN":
            logger.error(
                "Found a Sample on a supported platform with an unknown manufacturer.",
                sample=sample_object.id,
                platform_accession=sample_object.platform_accession_code,
                accession=sample_object.accession_code,
                manufacturer=sample_object.manufacturer,
                platform_name=sample_object.platform_name,
            )
            return ProcessorPipeline.NONE

    elif sample_object.has_raw:
        return ProcessorPipeline.SALMON

    # Shouldn't get here, but just in case
    return ProcessorPipeline.NONE


def determine_ram_amount(sample: Sample, job) -> int:
    """
    Determines the amount of RAM in MB required for a given ProcessorJob
    """

    if job.pipeline_applied == ProcessorPipeline.NO_OP.value:
        return 2048
    elif job.pipeline_applied == ProcessorPipeline.ILLUMINA_TO_PCL.value:
        return 2048
    elif job.pipeline_applied == ProcessorPipeline.AFFY_TO_PCL.value:
        platform = sample.platform_accession_code
        # Values via https://github.com/AlexsLemonade/refinebio/issues/54#issuecomment-373836510
        if "u133" in platform:
            return 2048
        if "gene" in platform:
            return 4096
        if "hta20" in platform or "huex10st" in platform:
            return 32768
        # Not sure what the ram usage of this platform is! Investigate!
        logger.debug("Unsure of RAM usage for platform! Using default.", platform=platform, job=job)
        return 2048
    elif job.pipeline_applied == ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL.value:
        return 2048
    elif job.pipeline_applied == ProcessorPipeline.AGILENT_ONECOLOR_TO_PCL.value:
        return 2048
    elif job.pipeline_applied == ProcessorPipeline.SALMON.value:
        return 12288
    elif job.pipeline_applied == ProcessorPipeline.NONE.value:
        return 1024
    else:
        logger.error(
            "Found a job without an expected pipeline!", job=job, pipeline=job.pipeline_applied
        )
        return 1024

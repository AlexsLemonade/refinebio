import random
import string

from enum import Enum, unique
from typing import List, Dict, Callable
from django.utils import timezone
from data_refinery_common.models import (
    ProcessorJob,
    Pipeline,
    Processor,
    Sample,
    OriginalFile,
    Dataset,
    ProcessorJobOriginalFileAssociation,
    ProcessorJobDatasetAssociation,
    OriginalFileSampleAssociation
)
from data_refinery_workers._version import __version__
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_worker_id, get_env_variable

logger = get_and_configure_logger(__name__)
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")


def start_job(job_context: Dict):
    """A processor function to start jobs.

    Record in the database that this job is being started and
    retrieves the job's batches from the database and adds them to the
    dictionary passed in with the key 'batches'.
    """
    job = job_context["job"]
    job.worker_id = get_worker_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    logger.info("Starting processor Job.", processor_job=job.id, pipeline=job.pipeline_applied)

    # Some jobs take OriginalFiles, other take Datasets
    if job.pipeline_applied not in ["SMASHER", "QN_REFERENCE"]:
        relations = ProcessorJobOriginalFileAssociation.objects.filter(processor_job=job)
        original_files = OriginalFile.objects.filter(id__in=relations.values('original_file_id'))

        if len(original_files) == 0:
            logger.error("No files found.", processor_job=job.id)
            job_context["success"] = False
            return job_context

        job_context["original_files"] = original_files
        original_file = job_context['original_files'][0]
        assocs = OriginalFileSampleAssociation.objects.filter(original_file=original_file)
        samples = Sample.objects.filter(id__in=assocs.values('sample_id')).distinct()
        job_context['samples'] = samples
        job_context["computed_files"] = []

    else:
        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)

        # This should never be more than one!
        dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
        dataset.is_processing = True
        dataset.save()

        # Get the samples to smash
        job_context["dataset"] = dataset
        job_context["samples"] = dataset.get_aggregated_samples()
        job_context["experiments"] = dataset.get_experiments()

        # Just in case
        job_context["original_files"] = []
        job_context["computed_files"] = []

    return job_context


def end_job(job_context: Dict, abort=False):
    """A processor function to end jobs.

    Record in the database that this job has completed and that
    the samples have been processed if not aborted.
    """
    job = job_context["job"]

    if "success" in job_context:
        success = job_context["success"]
    else:
        success = True

    if not abort:
        if job_context.get("success", False) and not (job_context["job"].pipeline_applied in ["SMASHER", "QN_REFERENCE"]):
            # This handles most of our cases
            for sample in job_context["samples"]:
                sample.is_processed = True
                sample.save()

            # Explicitly for the single-salmon scenario
            if 'sample' in job_context:
                sample = job_context['sample']
                sample.is_processed = True
                sample.save()

    # S3-sync Original Files
    for original_files in job_context['original_files']:
        # Ensure even distribution across S3 servers
        nonce = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
        result = original_files.sync_to_s3(S3_BUCKET_NAME, nonce + "_" + original_files.filename)
        if result:
            original_files.delete_local_file()

    # S3-sync Computed Files
    for computed_file in job_context['computed_files']:
        # Ensure even distribution across S3 servers
        nonce = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
        result = computed_file.sync_to_s3(S3_BUCKET_NAME, nonce + "_" + computed_file.filename)
        if result:
            computed_file.delete_local_file()

    # If the pipeline includes any steps, save it.
    pipeline = job_context['pipeline']
    if len(pipeline.steps):
        pipeline.save()

    job.success = success
    job.end_time = timezone.now()
    job.save()

    if success:
        logger.info("Processor job completed successfully.", processor_job=job.id)
    else:
        logger.info("Processor job failed!", processor_job=job.id)

    # Return Final Job context so testers can check it
    return job_context


def run_pipeline(start_value: Dict, pipeline: List[Callable]):
    """Runs a pipeline of processor functions.

    start_value must contain a key 'job_id' which is a valid id for a
    ProcessorJob record.

    Each processor fuction must accept a dictionary and return a
    dictionary.

    Any processor function which returns a dictionary containing a key
    of 'success' with a value of False will cause the pipeline to
    terminate with a call to utils.end_job.

    The key 'job' is reserved for the ProcessorJob currently being
    run.  The key 'batches' is reserved for the Batches that are
    currently being processed.  It is required that the dictionary
    returned by each processor function preserve the mappings for
    'job' and 'batches' that were passed into it.
    """
    job_id = start_value["job_id"]
    try:
        job = ProcessorJob.objects.get(id=job_id)
    except ProcessorJob.DoesNotExist:
        logger.error("Cannot find processor job record.", processor_job=job_id)
        return

    if len(pipeline) == 0:
        logger.error("Empty pipeline specified.",
                     procesor_job=job_id)

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        try:
            last_result = processor(last_result)
        except Exception:
            logger.exception("Unhandled exception caught while running processor function %s in pipeline",
                             processor.__name__,
                             processor_job=job_id)
            last_result["success"] = False
            return end_job(last_result)

        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor function %s failed. Terminating pipeline.",
                         processor.__name__,
                         processor_job=job_id,
                         failure_reason=last_result["job"].failure_reason)
            return end_job(last_result)

        if last_result.get("abort", False):
            return end_job(last_result, abort=True)

    return last_result


@unique
class PipelineEnum(Enum):
    """Hardcoded pipeline names."""

    AGILENT_TWOCOLOR = "Agilent Two Color"
    ARRAY_EXPRESS = "Array Express"
    ILLUMINA = "Illumina"
    NO_OP = "No Op"
    SALMON = "Salmon"
    SMASHER = "Smasher"
    TX_INDEX = "Transcriptome Index"
    QN_REFERENCE = "Quantile Normalization Reference"


@unique
class ProcessorEnum(Enum):
    """Hardcoded processor names in each pipeline."""

    # One processor in "Agilent Two Color" pipeline
    AGILENT_TWOCOLOR = "Agilent SCAN TwoColor"

    # One processor in "Array Express" pipeline
    AFFYMETRIX_SCAN = "Affymetrix SCAN"

    # One processor in "Illumina" pipeline
    ILLUMINA_SCAN = "Illumina SCAN"

    # One processor in "No Op" pipeline
    SUBMITTER_PROCESSED = "Submitter-processed"

    # Four processors in "Salmon" pipeline
    TXIMPORT = "Tximport"
    SALMON_QUANT = "Salmon Quant"
    MULTIQC = "MultiQC"
    SALMONTOOLS = "Salmontools"

    # No processors in "Smasher" pipeline (yet)

    # One processor in "Transcriptome Index" pipeline
    TX_INDEX = "Transcriptome Index"

    QN_REFERENCE = "Quantile Normalization Reference"


def createTestProcessors():
    """Creates dummy processors for all unit test cases.
    (This function should be called ONLY by test modules).
    """

    for label in ProcessorEnum:
        Processor.objects.create(name=label.value, version=__version__)

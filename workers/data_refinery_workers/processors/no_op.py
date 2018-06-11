import os
import shutil
import boto3
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputationalResult, ComputedFile, SampleResultAssociation
from data_refinery_common.utils import get_env_variable
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _no_op_processor_fn(job_context: Dict) -> Dict:
    """A processor which does nothing other than move files.

    Simply moves the file from its raw location to its
    processed location. Useful for handling data that has already been
    processed.
    """
    original_file = job_context["original_files"][0]

    # Create the output directory and path
    job_context["input_file_path"] = original_file.absolute_file_path
    base_directory, file_name = original_file.absolute_file_path.rsplit('/', 1)
    os.makedirs(base_directory + '/processed/', exist_ok=True)
    job_context["output_file_path"] = base_directory + '/processed/' + file_name

    # Copy the file to the new directory
    shutil.copyfile(job_context["input_file_path"], job_context["output_file_path"])

    # This is a NO-OP, but we make a ComputationalResult regardless.
    result = ComputationalResult()
    result.command_executed = "" # No op!
    result.is_ccdl = True
    result.system_version = __version__
    result.pipeline = "Submitter-processed"
    result.save()

    # Create a ComputedFile for the original file,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = job_context["output_file_path"]
        computed_file.filename = original_file.filename
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        # computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
        # TODO here: delete local file after S3 sync
        computed_file.save()
    except Exception:
        logger.error("Exception caught while moving file %s",
                         raw_path,
                         processor_job=job_context["job_id"])

        failure_reason = "Exception caught while moving file {}".format(file.name)
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context        

    for sample in job_context['samples']:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

    logger.info("Created %s", result)
    job_context["success"] = True
    return job_context


def no_op_processor(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _no_op_processor_fn,
                        utils.end_job])

import os
import shutil
import boto3
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputationalResult, ComputedFile
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

    sample = job_context["samples"][0]

    # This is a NO-OP, but we make a ComputationalResult regardless.
    result = ComputationalResult()
    result.command_executed = "" # No op!
    result.is_ccdl = True
    result.system_version = __version__
    result.save()

    # Create a ComputedFile for the sample,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = sample.source_absolute_file_path
        computed_file.filename = sample.source_filename
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
        # TODO here: delete local file after S3 sync
        computed_file.save()
    except Exception:
        logger.exception("Exception caught while moving file %s",
                         raw_path,
                         processor_job=job_context["job_id"])

        failure_reason = "Exception caught while moving file {}".format(file.name)
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context        

    logger.info("Created %s", result)
    job_context["success"] = True
    return job_context


def no_op_processor(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _no_op_processor_fn,
                        utils.end_job])

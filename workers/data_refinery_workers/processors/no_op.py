from typing import Dict
import os
import shutil
import boto3
from celery import shared_task
from data_refinery_common.models import batches, File
from data_refinery_workers.processors import utils

# import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def _no_op_processor_fn(job_context: Dict) -> Dict:
    """A processor which does nothing other than move files.

    Simply moves the batch's file from its raw location to its
    processed location. Useful for handling data that has already been
    processed.
    """
    file = File.objects.get(batch=job_context["batches"][0])
    raw_path = file.get_raw_path()

    try:
        if batches.USE_S3:
            client = boto3.client("s3")
            client.copy_object(Bucket=batches.S3_BUCKET_NAME,
                               Key=file.get_processed_path(),
                               CopySource={"Bucket": batches.S3_BUCKET_NAME,
                                           "Key": raw_path})
        else:
            os.makedirs(file.get_processed_dir(), exist_ok=True)
            shutil.copyfile(file.get_raw_path(),
                            file.get_processed_path())
    except Exception:
        logging.exception(("Exception caught while moving file %s for batch %d"
                           " during Job #%d."),
                          raw_path,
                          file.batch.id,
                          job_context["job_id"])
        failure_reason = "Exception caught while moving file {}".format(file.name)
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    try:
        file.remove_raw_files()
    except:
        # If we fail to remove the raw files, the job is still done
        # enough to call a success. However logging will be important
        # so the problem can be identified and the raw files cleaned up.
        logging.exception(("Exception caught while removing raw file %s for batch %d"
                           " during Job #%d."),
                          raw_path,
                          file.batch.id,
                          job_context["job_id"])

    job_context["success"] = True
    return job_context


@shared_task
def no_op_processor(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _no_op_processor_fn,
                        utils.end_job])

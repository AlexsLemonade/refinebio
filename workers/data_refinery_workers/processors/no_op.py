from typing import Dict
import shutil
import boto3
from celery import shared_task
from data_refinery_workers.processors import utils
from data_refinery_common import file_management
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


def no_op_processor_fn(job_context: Dict) -> Dict:
    """A processor which does nothing other than move files.

    Simply moves the batch's file from its raw location to its
    processed location. Useful for handling data that has already been
    processed.
    """
    batch = job_context["batches"][0]
    raw_path = file_management.get_raw_path(batch)

    try:
        if file_management.USE_S3:
            client = boto3.client("s3")
            client.copy_object(Bucket=file_management.S3_BUCKET_NAME,
                               Key=file_management.get_processed_path(batch),
                               CopySource={"Bucket": file_management.S3_BUCKET_NAME,
                                           "Key": raw_path})
        else:
            shutil.copyfile(file_management.get_raw_path(batch),
                            file_management.get_processed_path(batch))
    except Exception:
        logger.exception("Exception caught while moving file %s",
                         raw_path,
                         processor_job=job_context["job_id"],
                         batch=batch.id)

        failure_reason = "Exception caught while moving file {}".format(batch.name)
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    try:
        file_management.remove_raw_files(batch)
    except:
        # If we fail to remove the raw files, the job is still done
        # enough to call a success. However logging will be important
        # so the problem can be identified and the raw files cleaned up.
        logger.exception("Exception caught while removing raw file %s",
                         raw_path,
                         batch=batch.id,
                         processor_job=job_context["job_id"])

    job_context["success"] = True
    return job_context


@shared_task
def no_op_processor(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        no_op_processor_fn,
                        utils.end_job])

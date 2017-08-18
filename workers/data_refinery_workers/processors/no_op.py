from typing import Dict
import shutil
import boto3
from celery import shared_task
from data_refinery_workers.processors import utils
from data_refinery_common import file_management

# import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def no_op_processor_fn(kwargs: Dict) -> Dict:
    """A processor which does nothing other than move files.

    Simply moves the batch's file from its raw location to its
    processed location. Useful for handling data that has already been
    processed.
    """
    batch = kwargs["batches"][0]
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
        logging.exception(("Exception caught while moving file %s for batch %d"
                           " during Job #%d."),
                          raw_path,
                          batch.id,
                          kwargs["job_id"])
        kwargs["success"] = False
        return kwargs

    try:
        file_management.remove_raw_files(batch)
    except:
        # If we fail to remove the raw files, the job is still done
        # enough to call a success. However logging will be important
        # so the problem can be identified and the raw files cleaned up.
        logging.exception(("Exception caught while removing raw file %s for batch %d"
                           " during Job #%d."),
                          raw_path,
                          batch.id,
                          kwargs["job_id"])

    kwargs["success"] = True
    return kwargs


@shared_task
def no_op_processor(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        no_op_processor_fn,
                        utils.end_job])

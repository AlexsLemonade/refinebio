from __future__ import absolute_import, unicode_literals
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
from data_refinery_common import file_management
import logging

logger = get_task_logger(__name__)


def cel_to_pcl(kwargs: Dict) -> Dict:
    """Process .CEL files to .PCL format using R.

    Moves the .CEL file from the raw directory to the temp directory,
    calls ProcessCelFiles, uploads the processed file, and cleans up
    the raw/temp files. Because it uses the file_management module
    this works seamlessly whether S3 is being used or not.
    """
    # Array Express processor jobs have one batch per job.
    batch = kwargs["batches"][0]

    try:
        file_management.download_raw_file(batch)
    except Exception:
        logging.exception("Exception caught while retrieving %s for batch %d during Job #%d.",
                          file_management.get_raw_path(batch),
                          batch.id,
                          kwargs["job_id"])
        kwargs["success"] = False
        return kwargs

    temp_dir = file_management.get_temp_dir(batch)
    output_file = file_management.get_temp_post_path(batch)

    ro.r('source("/home/user/r_processors/cel_to_pcl.R")')
    ro.r['ProcessCelFiles'](
        temp_dir,
        "Hs",  # temporary until organism handling is more defined
        output_file)

    try:
        file_management.upload_processed_file(batch)
    except Exception:
        logging.exception(("Exception caught while uploading processed file %s for batch %d"
                           " during Job #%d."),
                          output_file,
                          batch.id,
                          kwargs["job_id"])
        kwargs["success"] = False
        return kwargs
    finally:
        file_management.remove_temp_directory(batch)

    try:
        file_management.remove_raw_files(batch)
    except:
        # If we fail to remove the raw files, the job is still done
        # enough to call a success. However logging will be important
        # so the problem can be identified and the raw files cleaned up.
        logging.exception(("Exception caught while removing raw files %s for batch %d"
                           " during Job #%d."),
                          output_file,
                          batch.id,
                          kwargs["job_id"])

    kwargs["success"] = True
    return kwargs


@shared_task
def affy_to_pcl(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        cel_to_pcl,
                        utils.end_job])

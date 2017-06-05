"""This processor is currently out of date. It is designed to process
multiple CEL files at a time, but that is not how we are going to
process Array Express files. I have rewritten the Array Express
surveyor/downloader to support this, but we don't have the new
processor yet. This will run, which is good enough for testing
the system, however since it will change so much the processor
itself is not yet tested.
"""


from __future__ import absolute_import, unicode_literals
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
from data_refinery_common import file_management
import logging

logger = get_task_logger(__name__)


def cel_to_pcl(kwargs: Dict):
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

    raw_file = file_management.get_temp_pre_path(batch)
    processed_file = file_management.get_temp_post_path(batch)

    # It's necessary to load the foreach library before calling SCANfast
    # because it doesn't load the library before calling functions
    # from it.
    ro.r("library('foreach')")

    ro.r['::']('SCAN.UPC', 'SCANfast')(
        raw_file,
        processed_file
    )

    try:
        file_management.upload_processed_file(batch)
    except Exception:
        logging.exception(("Exception caught while uploading processed file %s for batch %d"
                           " during Job #%d."),
                          processed_file,
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
        logging.exception(("Exception caught while uploading processed file %s for batch %d"
                           " during Job #%d."),
                          processed_file,
                          batch.id,
                          kwargs["job_id"])

    kwargs["success"] = True
    return kwargs


@shared_task
def affy_to_pcl(job_id):
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        cel_to_pcl,
                        utils.end_job])

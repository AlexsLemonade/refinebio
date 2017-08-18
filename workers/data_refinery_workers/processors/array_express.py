from __future__ import absolute_import, unicode_literals
import string
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
from data_refinery_common import file_management
import logging

logger = get_task_logger(__name__)


PACKAGE_NAME_CORRECTIONS = {
    "hugene10stv1hsentrezgprobe": "hugene10sthsentrezgprobe"
}


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

    input_file = file_management.get_temp_pre_path(batch)
    output_file = file_management.get_temp_post_path(batch)

    header = ro.r['::']('affyio', 'read.celfile.header')(input_file)

    # header is a list of vectors. [0][0] contains the package name.
    punctuation_table = str.maketrans(dict.fromkeys(string.punctuation))
    package_name = header[0][0].translate(punctuation_table).lower()
    # Headers can contain the version "v1" or "v2", which doesn't
    # appear in the brainarray package name. This replacement is
    # brittle, but the list of brainarray packages is relatively short
    # and we can monitor what packages are added to it and modify
    # accordingly. So far "v1" and "v2" are the only known versions
    # which must be accomodated in this way.
    package_name_without_version = package_name.replace("v1", "").replace("v2", "")
    brainarray_package = package_name_without_version + "hsentrezgprobe"

    # Prevents:
    # RRuntimeWarning: There were 50 or more warnings (use warnings()
    # to see the first 50)
    ro.r("options(warn=1)")

    # It's necessary to load the foreach library before calling SCANfast
    # because it doesn't load the library before calling functions
    # from it.
    ro.r("library('foreach')")

    ro.r['::']('SCAN.UPC', 'SCANfast')(
        input_file,
        output_file,
        probeSummaryPackage=brainarray_package
    )

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

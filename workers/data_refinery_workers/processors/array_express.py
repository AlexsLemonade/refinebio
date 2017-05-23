"""This processor is currently out of date. It is designed to process
multiple CEL files at a time, but that is not how we are going to
process Array Express files. I have rewritten the Array Express
surveyor/downloader to support this, but we don't have the new
processor yet. This will run, which is good enough for testing
the system, however since it will change so much the processor
itself is not yet tested."""


from __future__ import absolute_import, unicode_literals
import os
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils

logger = get_task_logger(__name__)


def cel_to_pcl(kwargs: Dict):
    # Array Express processor jobs have one batch per job.
    batch = kwargs["batches"][0]

    # TODO: the batch's name currently matches the file.
    # It should not store the file extension.
    raw_file = os.path.join(utils.ROOT_URI, "raw", batch.internal_location, batch.name)
    # raw_file = raw_file + "." + batch.raw_format

    target_directory = os.path.join(utils.ROOT_URI, "processed", batch.internal_location)
    os.makedirs(target_directory, exist_ok=True)

    processed_file = os.path.join(target_directory, batch.name)
    processed_file = processed_file + "." + batch.processed_format

    # It's necessary to load the foreach library before calling SCANfast
    # because it doesn't load the library before calling functions
    # from it.
    ro.r("library('foreach')")

    # rpy2 doesn't seem to support the double colon operator. Issue open here:
    # https://bitbucket.org/rpy2/rpy2/issues/408/accessing-a-librarys-function-via-the
    ro.r("library('SCAN.UPC')")
    ro.r["SCANfast"](
        raw_file,
        processed_file,
        probeSummaryPackage="hgu133plus2hsentrezgprobe"
    )

    return kwargs


@shared_task
def affy_to_pcl(job_id):
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        cel_to_pcl,
                        utils.end_job])

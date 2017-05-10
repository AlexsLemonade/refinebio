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

    from_directory = utils.ROOT_URI + "raw/" + batch.internal_location
    target_directory = utils.ROOT_URI + "processed/" + batch.internal_location
    os.makedirs(target_directory, exist_ok=True)
    new_name = batch.name + "." + batch.processed_format

    ro.r('source("/home/user/r_processors/process_cel_to_pcl.R")')
    ro.r['ProcessCelFiles'](
        from_directory,
        "Hs",  # temporary until organism handling is more defined
        target_directory + new_name)

    return kwargs


@shared_task
def affy_to_pcl(job_id):
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        cel_to_pcl,
                        utils.end_job])

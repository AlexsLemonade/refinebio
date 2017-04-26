from __future__ import absolute_import, unicode_literals
import os
import zipfile
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils

logger = get_task_logger(__name__)


def cel_to_pcl(kwargs: Dict):
    batch = kwargs["batch"]

    temp_directory = utils.ROOT_URI + "temp/" + batch.internal_location
    target_directory = utils.ROOT_URI + "processed/" + batch.internal_location
    os.makedirs(temp_directory, exist_ok=True)
    os.makedirs(target_directory, exist_ok=True)

    raw_file_name = batch.download_url.split('/')[-1]
    zip_location = (utils.ROOT_URI + "raw/" + batch.internal_location
                    + raw_file_name)

    # I need to figure out location tracking
    zip_ref = zipfile.ZipFile(zip_location, 'r')
    zip_ref.extractall(temp_directory)
    zip_ref.close()

    # Experiment code should be added to the batches data model
    experiment_code = raw_file_name.split('/')[0]
    new_name = experiment_code + ".pcl"

    ro.r('source("/home/user/r_processors/cel_to_pcl.r")')
    ro.r['ProcessCelFiles'](
        temp_directory,
        "Hs",  # temporary until organism discovery is working
        target_directory + new_name)

    return kwargs


@shared_task
def process_array_express(job_id):
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        cel_to_pcl,
                        # utils.cleanup_temp_data,
                        utils.end_job])

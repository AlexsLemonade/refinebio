from __future__ import absolute_import, unicode_literals
from typing import Dict
import rpy2.robjects as ro
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
# import rpy2
# print(rpy2.__version__)
logger = get_task_logger(__name__)


def pcl_to_cel(kwargs: Dict):
    batch = kwargs["batch"]

    ro.r('source("http://bioconductor.org/biocLite.R")')
    ro.r('library("affy")')
    ro.r('library("affyio")')
    ro.r('source("/home/user/process_to_PCL_brainarray.r")')
    ro.r('processCelFiles("/home/user/data_store/raw/ARRAY_EXPRESS/A-AFFY-1/E-MTAB-3050.raw.1.zip", 9606, /home/user/data_store/processed/ARRAY_EXPRESS/A-AFFY-1/E-MTAB-3050.pcl")')

    return kwargs


@shared_task
def process_array_express(job_id):
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        pcl_to_cel,
                        utils.end_job])

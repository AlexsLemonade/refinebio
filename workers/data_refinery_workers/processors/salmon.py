from __future__ import absolute_import, unicode_literals
# import string
from typing import Dict
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
# from data_refinery_common import file_management
# import logging
import subprocess


logger = get_task_logger(__name__)


def run_salmon(job_context: Dict) -> Dict:
    salmon = "salmon"
    index = "/home/user/mouse_index"
    read_1 = "/home/user/data_store/raw/IlluminaHiSeq2500/SALMON/ERR1680082_1.fastq"
    read_2 = "/home/user/data_store/raw/IlluminaHiSeq2500/SALMON/ERR1680082_2.fastq"
    output_dir = "/home/user/data_store/processed/IlluminaHiSeq2500/SALMON"
    cmd_str = "{} quant -l IU -i {} -1 {} -2 {} -p 20 -o {} --seqBias --gcBias --dumpEq"
    subprocess.Popen(cmd_str.format(salmon, index, read_1, read_2, output_dir).split())

    job_context["success"] = True
    return job_context


@shared_task
def salmon(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        run_salmon,
                        utils.end_job])

from __future__ import absolute_import, unicode_literals
import os
from typing import Dict
from celery import shared_task
from celery.utils.log import get_task_logger
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger
import subprocess


logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the file(s) from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    # Salmon processor jobs have only one batch per job, but may have
    # up to two files per batch.
    batch = job_context["batches"][0]
    files = batch.files

    for file in files:
        try:
            file.download_raw_file()
        except Exception:
            logger.exception("Exception caught while retrieving raw file %s",
                             file.get_raw_path(),
                             processor_job=job_context["job_id"],
                             batch=batch.id)

            failure_template = "Exception caught while retrieving raw file {}"
            job_context["job"].failure_reason = failure_template.format(file.name)
            job_context["success"] = False
            return job_context

    job_context["input_file_path"] = files[0].get_temp_pre_path()
    job_context["output_directory"] = os.path.join(files[0].get_temp_dir(), "output")

    if len(files) == 2:
        job_context["input_file_path_2"] = files[1].get_temp_pre_path()

    return job_context


def _zip_output(job_context: Dict) -> Dict:
    pass


def run_salmon(job_context: Dict) -> Dict:
    second_read_str = ""
    if "input_file_path_2" in job_context:
        second_read_str = " -2 {}".format(job_context["input_file_path_2"])

    # It appears that the index being missing didn't cause an error to
    # be thrown. I need to carefully test failures here.
    command_str = ("salmon quant -l IU -i {index} -1 {input_one_read_str} -p 20"
                   " -o {output_dir} --seqBias --gcBias --dumpEq --writeUnmappedNames")
    formatted_command = command_str.format(index="/home/user/data_store/mouse_index",
                                           input_one=job_context["input_file_path"],
                                           second_read_str=second_read_str,
                                           output_dir=job_context["output_directory"])
    subprocess.Popen(formatted_command.split())

    job_context["success"] = True
    return job_context


@shared_task
def salmon(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        run_salmon,
                        utils.end_job])

from __future__ import absolute_import, unicode_literals
import os
import tarfile
from typing import Dict
from celery import shared_task
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger
import subprocess


logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the file(s) from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_directory" to
    job_context so everything is prepared for processing. If the reads
    are paired then there will also be an "input_file_path_2" key
    added to job_context for the second read.
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
    # Salmon outputs an entire directory of files, so create a temp
    # directory to output it to until we can zip it to
    # files[0].get_temp_post_path()
    job_context["output_directory"] = os.path.join(files[0].get_temp_dir(), "output")
    os.makedirs(job_context["output_directory"], exist_ok=True)

    # Temporarily hardcoded until indices are implemented
    job_context["index_directory"] = "/home/user/data_store/mouse_index"

    if len(files) == 2:
        job_context["input_file_path_2"] = files[1].get_temp_pre_path()

    return job_context


def _zip_and_upload(job_context: Dict) -> Dict:
    file = job_context["batches"][0].files[0]
    # If there are paired reads... the file name will be based off of
    # the first file's name, but not the second.
    processed_path = file.get_temp_post_path()
    try:
        with tarfile.open(processed_path, "w:gz") as tar:
            tar.add(job_context["output_directory"], arcname=os.sep)
    except Exception:
        logger.exception("Exception caught while zipping processed directory %s",
                         job_context["output_directory"],
                         processor_job=job_context["job_id"],
                         batch=file.batch.id)

        file.remove_temp_directory()
        failure_template = "Exception caught while zipping processed directory {}"
        job_context["job"].failure_reason = failure_template.format(file.name)
        job_context["success"] = False
        return job_context

    try:
        file.upload_processed_file()
    except Exception:
        logger.exception("Exception caught while uploading processed file %s",
                         processed_path,
                         processor_job=job_context["job_id"],
                         batch=file.batch.id)

        file.remove_temp_directory()
        failure_template = "Exception caught while uploadiong processed file {}"
        job_context["job"].failure_reason = failure_template.format(processed_path)
        job_context["success"] = False
        return job_context

    file.remove_temp_directory()
    job_context["success"] = True
    return job_context


def _run_salmon(job_context: Dict) -> Dict:
    second_read_str = ""
    if "input_file_path_2" in job_context:
        second_read_str = " -2 {}".format(job_context["input_file_path_2"])

    command_str = ("salmon --no-version-check quant -l IU -i {index}"
                   " -1 {input_one}{second_read_str}"
                   " -p 20 -o {output_dir} --seqBias --gcBias --dumpEq --writeUnmappedNames")
    formatted_command = command_str.format(index=job_context["index_directory"],
                                           input_one=job_context["input_file_path"],
                                           second_read_str=second_read_str,
                                           output_dir=job_context["output_directory"])
    logger.info("Running Salmon Quant using the following shell command: %s",
                formatted_command,
                processor_job=job_context["job_id"],
                batch=job_context["batches"][0])

    completed_command = subprocess.run(formatted_command.split(),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)

    if completed_command.returncode == 1:
        stderr = str(completed_command.stderr)
        error_start = stderr.find("Error:")
        logger.error("Shell call to salmon failed with error message: %s",
                     stderr[error_start:],
                     processor_job=job_context["job_id"],
                     batch=job_context["batches"][0])

        job_context["batches"][0].files[0].remove_temp_directory()

        # The failure_reason column is only 256 characters wide.
        error_end = error_start + 200
        job_context["job"].failure_reason = ("Shell call to salmon failed because: "
                                             + stderr[error_start:error_end])
        job_context["success"] = False

    return job_context


@shared_task
def salmon(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _run_salmon,
                        _zip_and_upload,
                        utils.end_job])

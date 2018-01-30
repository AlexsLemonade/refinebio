from __future__ import absolute_import, unicode_literals
import os
import gzip
import subprocess
import tarfile
from typing import Dict
from data_refinery_workers.processors import utils
from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.models import BatchStatuses, BatchKeyValue
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "processor_job_"


def _set_job_prefix(job_context: Dict) -> Dict:
    job_context["job_dir_prefix"] = JOB_DIR_PREFIX + str(job_context["job_id"])
    return job_context


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
            file.download_raw_file(job_context["job_dir_prefix"])
        except Exception:
            logger.exception("Exception caught while retrieving raw file %s",
                             file.get_raw_path(),
                             processor_job=job_context["job_id"],
                             batch=batch.id)

            failure_template = "Exception caught while retrieving raw file {}"
            job_context["job"].failure_reason = failure_template.format(file.name)
            job_context["success"] = False
            return job_context

    num_files = len(files)
    if num_files > 2 or num_files < 1:
        failure_message = ("{} files were found for a Salmon job. There should never"
                           " be more than two.").format(str(num_files))
        logger.error(failure_message, batch=batch, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = failure_message
        job_context["success"] = False
        return job_context
    elif num_files == 2:
        job_context["input_file_path_2"] = files[1].get_temp_pre_path(
            job_context["job_dir_prefix"])

    job_context["input_file_path"] = files[0].get_temp_pre_path(job_context["job_dir_prefix"])
    # Salmon outputs an entire directory of files, so create a temp
    # directory to output it to until we can zip it to
    # files[0].get_temp_post_path()
    job_context["output_directory"] = os.path.join(
        files[0].get_temp_dir(job_context["job_dir_prefix"]), "output")
    os.makedirs(job_context["output_directory"], exist_ok=True)

    job_context["success"] = True
    return job_context


def _determine_index_length(job_context: Dict) -> Dict:
    """Determines whether to use the long or short salmon index.

    Adds the key 'kmer_size' to the job_context with a value of '23'
    if the short index is appropriate or '31' if the long index is
    appropriate. For more information on index length see the
    _create_index function of the transcriptome_index processor.
    """
    total_base_pairs = 0
    number_of_reads = 0
    counter = 1

    with gzip.open(job_context["input_file_path"], "rt") as input_file:
        for line in input_file:
            # In the FASTQ file format, there are 4 lines for each
            # read. Three of these contain metadata about the
            # read. The string representing the read itself is found
            # on the second line of each quartet.
            if counter % 4 == 2:
                total_base_pairs += len(line.replace("\n", ""))
                number_of_reads += 1

            counter += 1

    if "input_file_path_2" in job_context:
        with gzip.open(job_context["input_file_path_2"], "rt") as input_file:
            for line in input_file:
                if counter % 4 == 2:
                    total_base_pairs += len(line.replace("\n", ""))
                    number_of_reads += 1

                counter += 1

    if total_base_pairs / number_of_reads > 75:
        job_context["kmer_size"] = "31"
    else:
        job_context["kmer_size"] = "23"

    return job_context


def _download_index(job_context: Dict) -> Dict:
    """Downloads the appropriate Salmon Index for this batch.

    Salmon documentation states:

    "If you want to use Salmon in quasi-mapping-based mode, then you
    first have to build an Salmon index for your transcriptome."

    We have used the Data Refinery to build these indices already,
    this function retrieves the correct index for the organism and
    read length from Permanent Storage.
    """
    batch = job_context["batches"][0]
    try:
        index_batch = BatchKeyValue.objects.select_related("batch").filter(
            batch__source_type=Downloaders.TRANSCRIPTOME_INDEX.value,
            batch__organism_id=batch.organism_id,
            batch__status=BatchStatuses.PROCESSED.value,
            key="kmer_size",
            value=job_context["kmer_size"]
        ).all()[0].batch

        index_file = index_batch.files[0]
    except:
        logger.exception("Failed to find an index for organism %s with kmer_size of %s.",
                         batch.organism_name,
                         job_context["kmer_size"],
                         processor_job=job_context["job_id"],
                         batch=batch.id)

        failure_template = "Failed to find an index for organism {} with kmer_size of {}."
        job_context["job"].failure_reason = failure_template.format(batch.organism_name,
                                                                    job_context["kmer_size"])
        job_context["success"] = False
        return job_context

    temp_dir = batch.files[0].get_temp_dir(job_context["job_dir_prefix"])
    job_context["index_directory"] = os.path.join(temp_dir, "index")
    os.makedirs(job_context["index_directory"], exist_ok=True)
    try:
        # The index_batch file will be downloaded to a directory
        # specific to that batch.
        temp_path = index_file.download_processed_file(job_context["job_dir_prefix"])
        with tarfile.open(temp_path, "r:gz") as tarball:
            tarball.extractall(temp_dir)
    except:
        logger.exception("Failed to download and extract index tarball %s",
                         index_file.name,
                         index_batch=index_batch.id,
                         index_file=index_file.id,
                         processor_job=job_context["job_id"],
                         batch=batch.id)

        failure_template = "Failed to download and extract index tarball {}"
        job_context["job"].failure_reason = failure_template.format(index_file.name)
        job_context["success"] = False
        return job_context
    finally:
        # Remove temp dir for the index batch since we've extracted it to
        # the temp dir for this job.
        index_file.remove_temp_directory(job_context["job_dir_prefix"])

    job_context["success"] = True
    return job_context


def _zip_and_upload(job_context: Dict) -> Dict:
    """Zips the directory output by Salmon into a single file and uploads it.

    Adds the 'success' key to job_context because this function is the
    last in the job.
    """
    # If there are paired reads... the file name will be based off of
    # the first file's name, but not the second.
    # Note that this is a workaround for not having separate File
    # objects to represent processed files. Once that is fixed this
    # should simply use the processed file instead of one of the input
    # files.
    first_file = None
    for file in job_context["batches"][0].files:
        if file.name.find("_1."):
            first_file = file

    # If there is only a single read, just use that one.
    if first_file is None:
        first_file = job_context["batches"][0].files[0]

    processed_path = first_file.get_temp_post_path(job_context["job_dir_prefix"])
    try:
        with tarfile.open(processed_path, "w:gz") as tar:
            tar.add(job_context["output_directory"], arcname=os.sep)
    except Exception:
        logger.exception("Exception caught while zipping processed directory %s",
                         job_context["output_directory"],
                         processor_job=job_context["job_id"],
                         batch=first_file.batch.id)

        first_file.remove_temp_directory(job_context["job_dir_prefix"])
        failure_template = "Exception caught while zipping processed directory {}"
        job_context["job"].failure_reason = failure_template.format(first_file.name)
        job_context["success"] = False
        return job_context

    try:
        first_file.upload_processed_file(job_context["job_dir_prefix"])
    except Exception:
        logger.exception("Exception caught while uploading processed file %s",
                         processed_path,
                         processor_job=job_context["job_id"],
                         batch=first_file.batch.id)

        first_file.remove_temp_directory()
        failure_template = "Exception caught while uploading processed file {}"
        job_context["job"].failure_reason = failure_template.format(processed_path)
        job_context["success"] = False
        return job_context

    job_context["success"] = True
    return job_context


def _run_salmon(job_context: Dict) -> Dict:
    second_read_str = ""
    if "input_file_path_2" in job_context:
        second_read_str = " -2 {}".format(job_context["input_file_path_2"])

    command_str = ("salmon --no-version-check quant -l A -i {index}"
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

    ## To me, this looks broken: error codes are anything non-zero.
    ## However, Salmon (seems) to output with negative status codes
    ## even with succesful executions.
    ## Possibly related: https://github.com/COMBINE-lab/salmon/issues/55
    if completed_command.returncode == 1:
        stderr = str(completed_command.stderr)
        error_start = stderr.find("Error:")
        error_start = error_start if error_start != -1 else 0
        logger.error("Shell call to salmon failed with error message: %s",
                     stderr[error_start:],
                     processor_job=job_context["job_id"],
                     batch=job_context["batches"][0])

        job_context["batches"][0].files[0].remove_temp_directory(job_context["job_dir_prefix"])

        # The failure_reason column is only 256 characters wide.
        error_end = error_start + 200
        job_context["job"].failure_reason = ("Shell call to salmon failed because: "
                                             + stderr[error_start:error_end])
        job_context["success"] = False
    else:
        job_context["success"] = True

    return job_context


def salmon(job_id: int) -> None:
    """Main processor function for the Salmon Processor.

    Runs salmon quant command line tool, specifying either a long or
    short read length.
    """
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _set_job_prefix,
                        _prepare_files,
                        _determine_index_length,
                        _download_index,
                        _run_salmon,
                        _zip_and_upload,
                        utils.cleanup_raw_files,
                        utils.end_job])

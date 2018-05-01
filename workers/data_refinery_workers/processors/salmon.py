from __future__ import absolute_import, unicode_literals

import io
import json
import gzip
import os
import subprocess
import tarfile

from django.utils import timezone
from typing import Dict

import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError

from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    OrganismIndex, 
    ComputationalResult, 
    ComputationalResultAnnotation, 
    ComputedFile)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils

logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "processor_job_"
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
SKIP_PROCESSED = get_env_variable("SKIP_PROCESSED", True)


def _set_job_prefix(job_context: Dict) -> Dict:
    """ Sets the `job_dir_prefix` value in the job context object."""
    job_context["job_dir_prefix"] = JOB_DIR_PREFIX + str(job_context["job_id"])
    return job_context


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the file(s) from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_directory" to
    job_context so everything is prepared for processing. If the reads
    are paired then there will also be an "input_file_path_2" key
    added to job_context for the second read.
    """
    logger.debug("Preparing files..")

    original_files = job_context["original_files"]
    job_context["input_file_path"] = original_files[0].absolute_file_path
    if len(original_files) == 2:
        job_context["input_file_path_2"] = original_files[1].absolute_file_path

    # Salmon outputs an entire directory of files, so create a temp
    # directory to output it to until we can zip it to

    pre_part = original_files[0].absolute_file_path.split('/')[:-1]
    job_context["output_directory"] = '/'.join(pre_part) + '/processed/'
    os.makedirs(job_context["output_directory"], exist_ok=True)
    job_context["input_directory"] = '/'.join(pre_part) + '/'
    job_context["qc_directory"] = '/'.join(pre_part) + '/qc/'
    os.makedirs(job_context["qc_directory"], exist_ok=True)

    timestamp = str(timezone.now().timestamp()).split('.')[0]
    job_context["output_archive"] = '/'.join(pre_part) + '/result-' + timestamp +  '.tar.gz'
    os.makedirs(job_context["output_directory"], exist_ok=True)

    job_context['organism'] = job_context['original_files'][0].samples.first().organism
    job_context["success"] = True
    return job_context


def _determine_index_length(job_context: Dict) -> Dict:
    """Determines whether to use the long or short salmon index.

    Adds the key 'kmer_size' to the job_context with a value of '23'
    if the short index is appropriate or '31' if the long index is
    appropriate. For more information on index length see the
    _create_index function of the transcriptome_index processor.
    """
    logger.debug("Determining index length..")
    total_base_pairs = 0
    number_of_reads = 0
    counter = 1

    ### TODO: This process is really slow!
    ### Related: https://github.com/AlexsLemonade/refinebio/issues/157

    ### I bet there is a faster way of doing this, maybe by
    ### shelling and using UNIX!
    ### Python is single-core gunziping this line by line,
    ### I think it'd be faster to gunzip with multicores
    ### and then count lines that way, ex:

    ### $ pigz file.fasq.gz;
    ### $ cat file.fastq | echo $((`wc -l`/4))
    
    ### We may also be able to use io.BufferedReader to improve gzip speed
    ### Ex: (But does this work for text, not binary, gzip?)
    ### buffered_input = io.BufferedReader(input_file)
    ### for line in buffered_input.readlines():
    ###     if counter % 4 == 2:
    ###         total_base_pairs += len(line.replace("\n", ""))
    ###         number_of_reads += 1
    ###     counter += 1

    input_file = gzip.open(job_context["input_file_path"], "rt")

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
        job_context["index_length"] = "long"
    else:
        job_context["index_length"] = "short"

    return job_context


def _download_index(job_context: Dict) -> Dict:
    """Downloads the appropriate Salmon Index for this experiment.

    Salmon documentation states:

    "If you want to use Salmon in quasi-mapping-based mode, then you
    first have to build an Salmon index for your transcriptome."

    We have used the Data Refinery to build these indices already,
    this function retrieves the correct index for the organism and
    read length from Permanent Storage.
    """
    logger.debug("Downloading and installing index..")

    index_type = "TRANSCRIPTOME_" + job_context["index_length"].upper()
    index_object = OrganismIndex.objects.filter(organism=job_context['organism'], 
        index_type=index_type).order_by('created_at')[0]
    result = index_object.result
    files = ComputedFile.objects.filter(result=result)
    job_context["index_unpacked"] = '/'.join(files[0].absolute_file_path.split('/')[:-1])
    job_context["index_directory"] = job_context["index_unpacked"] + "/index"
    
    if not os.path.exists(job_context["index_directory"] + '/versionInfo.json'):
        with tarfile.open(files[0].absolute_file_path, "r:gz") as tarball:
            tarball.extractall(job_context["index_unpacked"])
    else:
        logger.info("Index already installed")

    job_context["success"] = True
    return job_context


def _zip_and_upload(job_context: Dict) -> Dict:
    """Zips the directory output by Salmon into a single file and uploads it.

    Adds the 'success' key to job_context because this function is the
    last in the job.
    """
    try:
        with tarfile.open(job_context['output_archive'], "w:gz") as tar:
            tar.add(job_context["output_directory"], arcname=os.sep)
    except Exception:
        logger.exception("Exception caught while zipping processed directory %s",
                         job_context["output_directory"],
                         processor_job=job_context["job_id"]
                        )
        failure_template = "Exception caught while zipping processed directory {}"
        job_context["job"].failure_reason = failure_template.format(first_file.name)
        job_context["success"] = False
        return job_context

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context["output_archive"]
    computed_file.filename = os.path.split(job_context["output_archive"])[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.is_public = True
    computed_file.result = job_context['result']
    computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
    # TODO here: delete local file after S3 sync#
    computed_file.save()

    job_context["success"] = True
    return job_context


def _run_salmon(job_context: Dict, skip_processed=SKIP_PROCESSED) -> Dict:
    """ """
    logger.debug("Running Salmon..")
    
    skip = False
    if skip_processed and os.path.exists(os.path.join(job_context['output_directory'] + 'quant.sf')):
        logger.info("Skipping pre-processed Salmon run!")
        skip = True

    # Salmon needs to be run differently for different sample types.
    # XXX: TODO: We need to tune the -p/--numThreads to the machines this process wil run on.
    if "input_file_path_2" in job_context:
        second_read_str = " -2 {}".format(job_context["input_file_path_2"])
        command_str = ("salmon --no-version-check quant -l A -i {index}"
                       " -1 {input_one}{second_read_str}"
                       " -p 20 -o {output_directory} --seqBias --gcBias --dumpEq --writeUnmappedNames")
        formatted_command = command_str.format(index=job_context["index_directory"],
                    input_one=job_context["input_file_path"],
                    second_read_str=second_read_str,
                    output_directory=job_context["output_directory"])
    else:
        # Related: https://github.com/COMBINE-lab/salmon/issues/83
        command_str = ("salmon --no-version-check quant -l A -i {index}"
               " -r {input_one}"
               " -p 20 -o {output_directory} --seqBias --dumpEq --writeUnmappedNames")
        formatted_command = command_str.format(index=job_context["index_directory"],
                    input_one=job_context["input_file_path"],
                    output_directory=job_context["output_directory"])

    logger.info("Running Salmon Quant using the following shell command: %s",
                formatted_command,
                processor_job=job_context["job_id"])

    job_context['time_start'] = timezone.now()
    if not skip:
        completed_command = subprocess.run(formatted_command.split(),
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
    job_context['time_end'] = timezone.now()

    ## To me, this looks broken: error codes are anything non-zero.
    ## However, Salmon (seems) to output with negative status codes
    ## even with succesful executions.
    ## Possibly related: https://github.com/COMBINE-lab/salmon/issues/55
    if not skip and completed_command.returncode == 1:

        stderr = str(completed_command.stderr)
        error_start = stderr.find("Error:")
        error_start = error_start if error_start != -1 else 0
        logger.error("Shell call to salmon failed with error message: %s",
                     stderr[error_start:],
                     processor_job=job_context["job_id"])

        # The failure_reason column is only 256 characters wide.
        error_end = error_start + 200
        job_context["job"].failure_reason = ("Shell call to salmon failed because: "
                                             + stderr[error_start:error_end])
        job_context["success"] = False
    else:

        result = ComputationalResult()
        result.command_executed = formatted_command
        result.system_version = __version__
        result.time_start = job_context['time_start']
        result.time_end = job_context['time_end']
        result.program_version = subprocess.run(['salmon', '--version'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).stderr.decode("utf-8").strip()
        result.is_ccdl = True
        result.save()

        with open(os.path.join(job_context['output_directory'], 'lib_format_counts.json')) as lfc_file:
            format_count_data = json.load(lfc_file)
            kv = ComputationalResultAnnotation()
            kv.data = format_count_data
            kv.result = result
            kv.is_public = True
            kv.save()
        with open(os.path.join(job_context['output_directory'], 'aux_info', 'meta_info.json')) as mi_file:
            meta_info = json.load(mi_file)
            kv = ComputationalResultAnnotation()
            kv.data = meta_info
            kv.result = result
            kv.is_public = True
            kv.save()

        job_context["result"] = result
        job_context["success"] = True

    return job_context

def _run_fastqc(job_context: Dict) -> Dict:
    """ Runs the `FastQC` package to generate the QC report.

    """

    command_str = ("./FastQC/fastqc --outdir={qc_directory} {files}")
    files = ' '.join(file.absolute_file_path for file in job_context['original_files'])
    formatted_command = command_str.format(qc_directory=job_context["qc_directory"],
                files=files)

    logger.info("Running FastQC using the following shell command: %s",
                formatted_command,
                processor_job=job_context["job_id"])

    completed_command = subprocess.run(formatted_command.split(),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)

    if completed_command.returncode != 0:
        stderr = str(completed_command.stderr)
        logger.error("Shell call to FastQC failed with error message: %s",
                     stderr,
                     processor_job=job_context["job_id"])

        # The failure_reason column is only 256 characters wide.
        job_context["job"].failure_reason = stderr[0:255]
        job_context["success"] = False

    return job_context

def _run_tximport(job_context: Dict) -> Dict:
    """ Runs the `tximport` R package to find genes. 
    
    For an R code example of tximport usage, see: 
    https://github.com/jaclyn-taroni/ref-txome/blob/79e2f64ffe6a71c5103a150bd3159efb784cddeb/4-athaliana_tximport.R

    Related: https://github.com/data-refinery/data-refinery/issues/82

    """

    try:
        tx_out = ro.r['::']('tximport', 'tximport')(job_context['output_directory'] + 'quant.sf', type="salmon", txOut=True)
    except RRuntimeError as e:
        # XXX/TODO: This doesn't work yet.
        # But, eventually, it will look something like this. ^^
        return job_context

    return job_context

def _run_salmontools(job_context: Dict) -> Dict:
    """ Run Salmontools to extract unmapped genes. 

    Related: https://github.com/data-refinery/data-refinery/issues/83
    """
    # XXX: TODO
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
                        _run_fastqc,


                        _determine_index_length,
                        _download_index,

                        _run_salmon,
                        _run_salmontools,
                        _run_tximport,
                        _run_fastqc,
                        _zip_and_upload,
                        utils.end_job])

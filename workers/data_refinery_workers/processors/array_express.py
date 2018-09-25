from __future__ import absolute_import, unicode_literals
import string
import warnings
from typing import Dict
import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError
from data_refinery_common.models import File
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the CEL file from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    # Array Express processor jobs have only one batch per job and one
    # file per batch.
    batch = job_context["batches"][0]
    file = File.objects.get(batch=batch)

    try:
        file.download_raw_file()
    except Exception:
        logger.exception("Exception caught while retrieving raw file %s",
                         file.get_raw_path(),
                         processor_job=job_context["job_id"],
                         batch=file.batch.id)

        failure_template = "Exception caught while retrieving raw file {}"
        job_context["job"].failure_reason = failure_template.format(file.name)
        job_context["success"] = False
        return job_context

    job_context["input_file_path"] = file.get_temp_pre_path()
    job_context["output_file_path"] = file.get_temp_post_path()
    return job_context


def _determine_brainarray_package(job_context: Dict) -> Dict:
    """Determines the right brainarray package to use for the file.

    Expects job_context to contain the key 'input_file'. Adds the key
    'brainarray_package' to job_context.
    """
    input_file = job_context["input_file_path"]
    try:
        header = ro.r['::']('affyio', 'read.celfile.header')(input_file)
    except RRuntimeError as e:
        # Array Express processor jobs have only one batch per job.
        file = File.objects.get(batch=job_context["batches"][0])
        file.remove_temp_directory()

        error_template = ("Unable to read Affy header in input file {0}"
                          " while running AFFY_TO_PCL due to error: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.error(error_message, processor_job=job_context["job"].id)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        return job_context

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
    job_context["brainarray_package"] = package_name_without_version + "hsentrezgprobe"
    return job_context


def _run_scan_upc(job_context: Dict) -> Dict:
    """Processes an input CEL file to an output PCL file.

    Does so using the SCAN.UPC package's SCANfast method using R.
    Expects job_context to contain the keys 'input_file', 'output_file',
    and 'brainarray_package'.
    """
    input_file = job_context["input_file_path"]

    try:
        # It's necessary to load the foreach library before calling SCANfast
        # because it doesn't load the library before calling functions
        # from it.
        ro.r("suppressMessages(library('foreach'))")

        # Prevents:
        # RRuntimeWarning: There were 50 or more warnings (use warnings()
        # to see the first 50)
        ro.r("options(warn=1)")

        # All R messages are turned into Python 'warnings' by rpy2. By
        # filtering all of them we silence a lot of useless output
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scan_upc = ro.r['::']('SCAN.UPC', 'SCANfast')
            scan_upc(input_file,
                     job_context["output_file_path"],
                     probeSummaryPackage=job_context["brainarray_package"])

    except RRuntimeError as e:
        error_template = ("Encountered error in R code while running AFFY_TO_PCL"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        # Array Express processor jobs have only one batch per job and
        # one file per batch
        file = File.objects.get(batch=job_context["batches"][0])
        file.remove_temp_directory()

    return job_context


def affy_to_pcl(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _determine_brainarray_package,
                        _run_scan_upc,
                        utils.upload_processed_files,
                        utils.cleanup_raw_files,
                        utils.end_job])

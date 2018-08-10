from __future__ import absolute_import, unicode_literals

import os
import string
import warnings
from django.utils import timezone
from typing import Dict

import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    OriginalFile,
    ComputationalResult,
    ComputedFile,
    SampleResultAssociation,
    SampleComputedFileAssociation,
    Processor,
    Pipeline
)

from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.utils import get_env_variable
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the CEL file from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    original_file = job_context["original_files"][0]
    job_context["input_file_path"] = original_file.get_synced_file_path()
    # This is ugly, I'm sorry.
    # Turns /home/user/data_store/E-GEOD-8607/raw/foo.cel into /home/user/data_store/E-GEOD-8607/processed/foo.cel
    pre_part = original_file.absolute_file_path.split('/')[:-2]
    end_part = original_file.absolute_file_path.split('/')[-1]
    job_context["output_file_path"] = '/'.join(pre_part) + '/processed/' + end_part
    job_context["output_file_path"] = job_context["output_file_path"].replace('.CEL', '.PCL')
    job_context["output_file_path"] = job_context["output_file_path"].replace('.cel', '.PCL')

    return job_context


def _create_ensg_pkg_map() -> Dict:
    """Reads the text file that was generated when installing ensg R
    packages, and returns a map whose keys are chip names and values are
    the corresponding BrainArray ensg package name.
    """
    ensg_pkg_filename = "/home/user/r_ensg_probe_pkgs.txt"
    chip2pkg = dict()
    with open(ensg_pkg_filename) as file_handler:
        for line in file_handler:
            tokens = line.strip("\n").split("\t")
            # tokens[0] is (normalized) chip name,
            # tokens[1] is the package's URL in this format:
            # http://mbni.org/customcdf/<version>/ensg.download/<pkg>_22.0.0.tar.gz
            pkg_name = tokens[1].split("/")[-1].split("_")[0]
            chip2pkg[tokens[0]] = pkg_name

    return chip2pkg


def _determine_brainarray_package(job_context: Dict) -> Dict:
    """Determines the right brainarray package to use for the file.

    Expects job_context to contain the key 'input_file'. Adds the key
    'brainarray_package' to job_context.
    """
    input_file = job_context["input_file_path"]
    try:
        header = ro.r['::']('affyio', 'read.celfile.header')(input_file)
    except RRuntimeError as e:
        error_template = ("Unable to read Affy header in input file {0}"
                          " while running AFFY_TO_PCL due to error: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.error(error_message, processor_job=job_context["job"].id)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        return job_context

    # header is a list of vectors. [0][0] contains the package name.
    punctuation_table = str.maketrans(dict.fromkeys(string.punctuation))
    # Normalize header[0][0]
    package_name = header[0][0].translate(punctuation_table).lower()

    # Headers can contain the version "v1" or "v2", which doesn't
    # appear in the brainarray package name. This replacement is
    # brittle, but the list of brainarray packages is relatively short
    # and we can monitor what packages are added to it and modify
    # accordingly. So far "v1" and "v2" are the only known versions
    # which must be accomodated in this way.
    # Related: https://github.com/data-refinery/data-refinery/issues/141
    package_name_without_version = package_name.replace("v1", "").replace("v2", "")
    chip_pkg_map = _create_ensg_pkg_map()
    job_context["brainarray_package"] = chip_pkg_map.get(package_name_without_version, None)
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
            job_context['time_start'] = timezone.now()

            # Related: https://github.com/AlexsLemonade/refinebio/issues/64
            if job_context["brainarray_package"]:
                scan_upc(input_file,
                         job_context["output_file_path"],
                         probeSummaryPackage=job_context["brainarray_package"])
            else:
                scan_upc(input_file,
                         job_context["output_file_path"])
            job_context['time_end'] = timezone.now()

    except RRuntimeError as e:
        error_template = ("Encountered error in R code while running AFFY_TO_PCL"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

    return job_context

def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    result = ComputationalResult()
    result.commands.append('SCAN.UPC::SCANfast')
    result.is_ccdl = True
    result.is_public = True

    result.time_start = job_context['time_start']
    result.time_end = job_context['time_end']
    result.pipeline = "Affymetrix SCAN"  # TODO: should be removed
    try:
        processor_key = "AFFYMETRIX_SCAN"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return handle_processor_exception(job_context, processor_key, e)

    result.save()
    job_context['pipeline'].steps.append(result.id)

    # Create a ComputedFile for the sample,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = job_context["output_file_path"]
        computed_file.filename = os.path.split(job_context["output_file_path"])[-1]
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        computed_file.save()
        job_context['computed_files'].append(computed_file)
    except Exception:
        logger.exception("Exception caught while saving file %s to S3",
                         computed_file.filename,
                         processor_job=job_context["job_id"],
                         )
        failure_reason = "Exception caught while saving file to S3"
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    for sample in job_context['samples']:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample,
            computed_file=computed_file)

    logger.info("Created %s", result,
                processor_job=job_context["job_id"])
    job_context["success"] = True

    return job_context

def affy_to_pcl(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.ARRAY_EXPRESS.value)
    utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                       [utils.start_job,
                        _prepare_files,
                        _determine_brainarray_package,
                        _run_scan_upc,
                        _create_result_objects,
                        # utils.upload_processed_files,
                        # utils.cleanup_raw_files,
                        utils.end_job])

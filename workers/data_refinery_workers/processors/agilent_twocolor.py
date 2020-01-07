from __future__ import absolute_import, unicode_literals

import os
import warnings
from typing import Dict

from django.utils import timezone

import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    OriginalFile,
    Pipeline,
    Processor,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Populate our job_context with appropriate inputs and outputs

    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    original_file = job_context["original_files"][0]
    job_context["input_file_path"] = original_file.absolute_file_path
    # Turns /home/user/data_store/E-GEOD-8607/raw/foo.txt into /home/user/data_store/E-GEOD-8607/processed/foo.cel
    pre_part = original_file.absolute_file_path.split("/")[:-2]
    end_part = original_file.absolute_file_path.split("/")[-1]

    os.makedirs("/".join(pre_part) + "/processed/", exist_ok=True)
    job_context["output_file_path"] = "/".join(pre_part) + "/processed/" + end_part
    job_context["output_file_path"] = job_context["output_file_path"].replace(".txt", ".PCL")

    return job_context


def _run_scan_twocolor(job_context: Dict) -> Dict:
    """Processes an input TXT file to an output PCL file.

    Does so using the SCAN.UPC package's SCAN_TwoColor method using R.
    Expects job_context to contain the keys 'input_file' and 'output_file'.
    """
    input_file = job_context["input_file_path"]

    try:
        # It's necessary to load the foreach library before calling SCAN_TwoColor
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
            scan_upc = ro.r["::"]("SCAN.UPC", "SCAN_TwoColor")
            job_context["time_start"] = timezone.now()

            # XXX: Current bug lives here. See below.
            scan_upc(input_file, job_context["output_file_path"])
            job_context["time_end"] = timezone.now()

    except RRuntimeError as e:

        # There is a bug in Bioconductor that causes this, but it happens
        # _after_ our output file is made and written to disk so..
        # we can ignore it!
        # See it for yourself here: https://github.com/Miserlou/SCAN.UPC/blob/master/R/TwoColor.R#L122
        bugstring = 'Error in sampleNames(expressionSet) = sampleNames : \n  could not find function "sampleNames<-"\n'
        if str(e) == bugstring:
            job_context["time_end"] = timezone.now()
            job_context["success"] = True
        else:
            error_template = (
                "Encountered error in R code while running AGILENT_TWOCOLOR_TO_PCL"
                " pipeline during processing of {0}: {1}"
            )
            error_message = error_template.format(input_file, str(e))

            logger.error(error_message, processor_job=job_context["job_id"])
            job_context["job"].failure_reason = error_message
            job_context["success"] = False

    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    result = ComputationalResult()
    result.commands.append("SCAN.UPC::SCAN_TwoColor")
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    try:
        processor_key = "AGILENT_TWOCOLOR"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)

    result.save()
    job_context["pipeline"].steps.append(result.id)

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
        job_context["computed_files"].append(computed_file)
    except Exception:
        logger.exception(
            "Exception caught while moving file %s to S3",
            computed_file.filename,
            processor_job=job_context["job_id"],
        )
        failure_reason = "Exception caught while moving file to S3"
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    for sample in job_context["samples"]:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample, computed_file=computed_file
        )

    logger.info("Created %s", result)
    job_context["success"] = True

    return job_context


def agilent_twocolor_to_pcl(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.AGILENT_TWOCOLOR.value)
    utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _prepare_files,
            _run_scan_twocolor,
            _create_result_objects,
            utils.end_job,
        ],
    )

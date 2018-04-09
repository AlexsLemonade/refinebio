from __future__ import absolute_import, unicode_literals

import csv
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
    Sample
)
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.utils import get_env_variable
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the TXT file from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    original_file = job_context["original_files"][0]
    job_context["input_file_path"] = original_file.absolute_file_path
    # Turns /home/user/data_store/E-GEOD-8607/raw/foo.txt into /home/user/data_store/E-GEOD-8607/processed/foo.cel
    pre_part = original_file.absolute_file_path.split('/')[:-2]
    end_part = original_file.absolute_file_path.split('/')[-1]
    job_context["output_file_path"] = '/'.join(pre_part) + '/processed/' + end_part
    job_context["output_file_path"] = job_context["output_file_path"].replace('.txt', '.PCL')

    return job_context

def _collect_samples(job_context: Dict) -> Dict:
    """
    Detection of column type in the data depends upon knowing the sample names in advance.
    """

    original_files = job_context['original_files']
    samples = Sample.objects.filter(id__in=original_files.values('sample_id'))
    job_context['samples'] = samples

    return job_context

def _detect_columns(job_context: Dict) -> Dict:
    """ Detect which columns match to which inputs.

    Related: https://github.com/AlexsLemonade/refinebio/issues/86#issuecomment-379308817

    We need to find:

        First column should be ID_REF or PROBE_ID and the type should be string.
        Detection Pval column
        Expression column (contains sample title and NOT 'BEAD')

    """

    input_file = job_context["input_file_path"]
    headers = None
    with open(input_file, 'r') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            headers = row
            break

    import pdb
    pdb.set_trace()

    # First the probe ID column
    if headers[0] not in ['ID_REF', 'PROBE_ID']:
        job_context["job"].failure_reason = "Could not find ID reference column"
        job_context["success"] = False
        return job_context
    else:
        job_context['probeId'] = headers[0]


    return job_context

def _run_illumina(job_context: Dict) -> Dict:
    """Processes an input TXT file to an output PCL file.

    Does so using a custom R script.
    Expects job_context to contain the keys 'input_file', 'output_file'.
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
            scan_upc = ro.r['::']('XXX', 'XXX')
            job_context['time_start'] = timezone.now()
            scan_upc(input_file,
                     job_context["output_file_path"])
            job_context['time_end'] = timezone.now()

    except RRuntimeError as e:
        error_template = ("Encountered error in R code while running AGILENT_TWOCOLOR_TO_PCL"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(input_file, str(e))

        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

    return job_context

def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    result = ComputationalResult()
    result.command_executed = "Illumina.R" # Need a better way to represent this R code.
    result.is_ccdl = True
    result.is_public = True
    result.system_version = __version__
    result.program_version = "XXX"
    result.time_start = job_context['time_start']
    result.time_end = job_context['time_end']
    result.save()

    # Create a ComputedFile for the sample,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = job_context["output_file_path"]
        computed_file.filename = os.path.split(job_context["output_file_path"])[-1]
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        # computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
        # TODO here: delete local file after S3 sync
        computed_file.save()
    except Exception:
        logger.exception("Exception caught while moving file %s to S3",
                         computed_file.filename,
                         processor_job=job_context["job_id"],
                         )
        failure_reason = "Exception caught while moving file to S3"
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    logger.info("Created %s", result)
    job_context["success"] = True

    return job_context

def illumina_to_pcl(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _collect_samples,
                        _detect_columns,
                        _run_illumina,
                        _create_result_objects,
                        utils.end_job])

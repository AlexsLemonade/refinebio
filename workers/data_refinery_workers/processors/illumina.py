from __future__ import absolute_import, unicode_literals

import csv
import os
import string
import subprocess
import multiprocessing
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
    Sample,
    OriginalFileSampleAssociation,
    SampleResultAssociation
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

    # Sanitize this file so R doesn't choke.
    # Some have comments, some have non-comment-comments.
    file_input = open(job_context["input_file_path"], "r")
    file_output = open(job_context["input_file_path"] + ".sanitized", "w")
    for line in file_input:
        if '#' not in line and \
        line.strip() != '' and \
        line != '\n' and \
        '\t' in line and \
        line[0] != '\t' and \
        line != None:
            file_output.write(line)    
    file_input.close()
    file_output.close()
    job_context['input_file_path'] = job_context["input_file_path"] + ".sanitized"

    return job_context

def _detect_columns(job_context: Dict) -> Dict:
    """ Detect which columns match to which inputs.

    Related: https://github.com/AlexsLemonade/refinebio/issues/86#issuecomment-379308817

    We need to find:

        First column should be ID_REF or PROBE_ID and the type should be string.
        Detection Pval column
        Expression column (contains sample title and NOT 'BEAD')

        Header examples:     
            ['ID_REF', 'LV-C&si-Control-1', 'Detection Pval',
            'LV-C&si-Control-2', 'Detection Pval', 'LV-C&si-Control-3', 'Detection
            Pval', 'LV-C&si-EZH2-1', 'Detection Pval', 'LV-C&si-EZH2-2', 'Detection
            Pval', 'LV-C&si-EZH2-3', 'Detection Pval', 'LV-EZH2&si-EZH2-1',
            'Detection Pval', 'LV-EZH2&si-EZH2-2', 'Detection Pval', 'LV-EZH2&si-
            EZH2-3', 'Detection Pval', 'LV-T350A&si-EZH2-1', 'Detection Pval', 'LV-
            T350A&si-EZH2-2', 'Detection Pval', 'LV-T350A&si-EZH2-3', 'Detection
            Pval']

    """

    try:
        input_file = job_context["input_file_path"]
        headers = None
        with open(input_file, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            for row in tsvin:

                # Skip sparse header row
                if row[0] == "":
                    continue

                # Skip comment rows
                joined = ''.join(row)
                if '#' in joined or '' == joined:
                    continue

                headers = row
                break

        # First the probe ID column
        if headers[0] not in ['ID_REF', 'PROBE_ID']:
            job_context["job"].failure_reason = "Could not find ID reference column"
            job_context["success"] = False
            return job_context
        else:
            job_context['probeId'] = headers[0]

        # Then the detection Pvalue string, which is always(?) some form of 'Detection Pval'
        for header in headers:
            if 'DETECTION_PVAL' in header.upper().replace(' ', '_'):
                # Could be <SAMPLE_ID>.Detection Pval, etc
                if '.' in header:
                    job_context['detectionPval'] = header.split('.')[-1]
                else:
                    job_context['detectionPval'] = header
                break
        else:
            job_context["job"].failure_reason = "Could not detect PValue column!"
            job_context["success"] = False
            return job_context

        # Then, finally, create an absolutely bonkers regular expression
        # which will explictly hit on any sample which contains a 4sample
        # ID _and_ ignores the magical word 'BEAM'. Great!
        # TODO: What to do about the ones that say `.AVG_Signal` ?
        column_ids = ""
        for sample in job_context['samples']:
            for offset, header in enumerate(headers, start=1):

                if sample.title == header:
                    column_ids = column_ids + str(offset) + "," 
                    continue
                if header.upper().replace(' ', '_') == "RAW_VALUE":
                    column_ids = column_ids + str(offset) + "," 
                    continue
                if sample.title in header and \
                'BEAD' not in header.upper() and \
                'NARRAYS' not in header.upper() and \
                'ARRAY_STDEV' not in header.upper() and \
                'PVAL' not in header.upper().replace(' ', '').replace('_', ''):
                    column_ids = column_ids + str(offset) + "," 
                    continue
        for offset, header in enumerate(headers, start=1):
            if 'AVG_Signal' in header:
                column_ids = column_ids + str(offset) + "," 
                continue

        column_ids = column_ids[:-1]
        job_context['columnIds'] = column_ids
    except Exception as e:
        job_context["job"].failure_reason = str(e)
        job_context["success"] = False
        return job_context

    return job_context

def _run_illumina(job_context: Dict) -> Dict:
    """Processes an input TXT file to an output PCL file.

    Does so using a custom R script.
    Expects job_context to contain the keys 'input_file_path', 'output_file_path'.
    """
    input_file_path = job_context["input_file_path"]

    # We need to detect the sub-platform. Looking for the chip version
    # in the metadata is the best we can do for right now.
    annotation = job_context['samples'][0].sampleannotation_set.all()[0]
    annotation_data = str(annotation.data).encode('utf-8').upper()
    if "V2".encode() in annotation_data or "V 2".encode() in annotation_data:
        platform = "illuminaHumanv2"
    else:
        platform = "illuminaHumanv4"

    try:
        job_context['time_start'] = timezone.now()

        result = subprocess.check_output([
                "/usr/bin/Rscript", 
                "--vanilla", 
                "/home/user/data_refinery_workers/processors/illumina.R",
                "--probeId", job_context['probeId'],
                "--expression", job_context['columnIds'],
                "--detection", job_context['detectionPval'],
                "--platform", platform,
                "--inputFile", job_context['input_file_path'],
                "--outputFile", job_context['output_file_path'],
                "--cores", str(multiprocessing.cpu_count())
            ])

        job_context['time_end'] = timezone.now()

    except Exception as e:
        error_template = ("Encountered error in R code while running illumina.R"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(job_context['input_file_path'], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

    return job_context

def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    result = ComputationalResult()
    result.command_executed = "illumina.R" # Need a better way to represent this R code.
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

        logger.info("Created new ComputedFile: " + computed_file.absolute_file_path + " [" + str(computed_file.sha1) + "](" + str(computed_file.size_in_bytes) + ")",
            computed_file=computed_file)

    except Exception:
        logger.exception("Exception caught while moving file %s to S3",
                         computed_file.filename,
                         processor_job=job_context["job_id"],
                         )
        failure_reason = "Exception caught while moving file to S3"
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    for sample in job_context['samples']:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

    logger.info("Created %s", result)
    job_context["success"] = True

    return job_context

def illumina_to_pcl(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _detect_columns,
                        _run_illumina,
                        _create_result_objects,
                        utils.end_job])

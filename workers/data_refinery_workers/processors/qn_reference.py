from __future__ import absolute_import, unicode_literals

import os
import subprocess
import random
import string
import time
import warnings
from django.utils import timezone
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    OriginalFile,
    ComputationalResult,
    ComputedFile,
    SampleResultAssociation,
    SampleComputedFileAssociation,
    ComputationalResultAnnotation,
    Processor,
    Pipeline
)

from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils, smasher

from data_refinery_common.utils import get_env_variable
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)

def _prepare_input(job_context: Dict) -> Dict:
    
    # We're going to use the smasher outside of the smasher.
    # I'm not crazy about this yet. Maybe refactor later,
    # but I need the data now.

    job_context = smasher._prepare_files(job_context)
    job_context = smasher._smash(job_context)

    if not 'final_frame' in job_context.keys():
        logger.error("Unable to prepare files for creating QN target.", 
            job_id=job_context['job'].id)
        job_context['success'] = False
        return job_context

    # We only need the resulting frame, not the entire archive
    outfile = '/tmp/' + str(time.time()).split('.')[0] + '.tsv'
    job_context['final_frame'].to_csv(outfile, sep='\t', encoding='utf-8')
    job_context['smashed_file'] = outfile
    job_context['target_file'] = outfile + '.target.tsv'

    return job_context

def _quantile_normalize(job_context: Dict) -> Dict:

    try:
        job_context['time_start'] = timezone.now()

        result = subprocess.check_output([
                "/usr/bin/Rscript",
                "--vanilla",
                "/home/user/data_refinery_workers/processors/qn_reference.R",
                "--inputFile", job_context['smashed_file'],
                "--outputFile", job_context['target_file']
            ])

        job_context['time_end'] = timezone.now()

    except Exception as e:
        error_template = ("Encountered error in R code while running qn_reference.R"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(job_context['smashed_file'], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

    except RRuntimeError as e:
        error_template = ("Encountered error in R code while running qn_reference"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False


    return job_context

def _create_result_objects(job_context: Dict) -> Dict:

    result = ComputationalResult()
    result.commands.append("qn_reference.R") # Need a better way to represent this R code.
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context['time_start']
    result.time_end = job_context['time_end']
    try:
        processor_key = "QN_REFERENCE"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context['target_file']
    computed_file.filename = job_context['target_file'].split('/')[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.is_smashable = False
    computed_file.is_qn_target = True
    computed_file.result = result
    computed_file.save()

    annotation = ComputationalResultAnnotation()
    annotation.result = result
    annotation.data = {
        "organism_id": job_context['samples']['ALL'][0].organism_id,
        "platform_accession_code": job_context['samples']['ALL'][0].platform_accession_code
    }
    annotation.save()

    job_context['result'] = result
    job_context['computed_files'] = [computed_file]
    job_context['success'] = True
    return job_context

def create_qn_reference(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.QN_REFERENCE.value)
    job_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                       [utils.start_job,
                        _prepare_input,
                        _quantile_normalize,
                        _create_result_objects,
                        utils.end_job])
    return job_context

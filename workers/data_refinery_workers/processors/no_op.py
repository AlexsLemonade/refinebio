import csv
import os
import shutil
import boto3

import subprocess
import numpy as np
import pandas as pd
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputationalResult, ComputedFile, SampleResultAssociation
from data_refinery_common.utils import get_env_variable, get_internal_microarray_accession
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """A processor which does nothing other than move files.

    Simply moves the file from its raw location to its
    processed location. Useful for handling data that has already been
    processed.
    """
    original_file = job_context["original_files"][0]

    # Create the output directory and path
    job_context["input_file_path"] = original_file.absolute_file_path
    job_context["input_file_path"] = original_file.absolute_file_path
    base_directory, file_name = original_file.absolute_file_path.rsplit('/', 1)
    os.makedirs(base_directory + '/processed/', exist_ok=True)
    job_context["output_file_path"] = base_directory + '/processed/' + file_name

    # XXX: Are there files we still want to do this to?
    # Copy the file to the new directory
    # shutil.copyfile(job_context["input_file_path"], job_context["output_file_path"])
    # job_context["success"] = True

    # Make sure header column is correct
    with open(job_context["input_file_path"], 'r') as tsv_in:
        tsv_in = csv.reader(tsv_in, delimiter='\t')
        for line in tsv_in:
            row = line
            joined = ''.join(row)
            break

    # If ID_REF already, we're good.
    if 'ID_REF' not in joined:
        try:
            float(row[1])
            # Okay, there's no header so can just prepend to the file.
            with open(job_context["input_file_path"], 'r') as f:
                all_content = f.read()

            job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
            with open(job_context["input_file_path"], 'w+') as f:
                f.seek(0, 0)
                f.write('ID_REF\tVALUE' + '\n' + all_content)

        except ValueError:
            # There is already a header row. Let's replace it with ID_Ref
            with open(job_context["input_file_path"], 'r') as f:
                all_content = f.read()

            job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
            with open(job_context["input_file_path"], 'w+') as f:
                f.seek(0, 0)
                f.write('ID_REF\tVALUE' + '\n' + ("\n".join(all_content[1:])))
            job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"

    # Platform
    job_context["platform"] = job_context["samples"][0].platform_accession_code
    job_context["internal_accession"] = get_internal_microarray_accession(job_context["platform"])
    if not job_context["internal_accession"]:
        logger.error("Failed to find internal accession for code", code=job_context["platform"])
        job_context['success'] = False

    return job_context


def _convert_genes(job_context: Dict) -> Dict:
    """ Convert to Ensembl genes if we can"""

    gene_index_path = "/home/user/data_store/" + job_context["internal_accession"] + ".tsv.gz"
    if not os.path.exists(gene_index_path):
        logger.error("Missing gene index file for platform!". ,
            platform=job_context["internal_accession"])
        job_context["success"] = False
        return job_context

    try:
        result = subprocess.check_output([
                "/usr/bin/Rscript", 
                "--vanilla", 
                "/home/user/data_refinery_workers/processors/gene_convert.R",
                "--platform", job_context["internal_accession"],
                "--inputFile", job_context['input_file_path'],
                "--outputFile", job_context['output_file_path']
            ])
    except Exception as e:
        error_template = ("Encountered error in R code while running gene_convert.R"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(job_context['input_file_path'], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        return job_context

    job_context["success"] = True
    return job_context


def _create_result(job_context: Dict) -> Dict:
    """ Create the actual Result object"""

    # This is a NO-OP, but we make a ComputationalResult regardless.
    result = ComputationalResult()
    result.command_executed = "gene_convert.R"
    result.is_ccdl = True
    result.system_version = __version__
    result.pipeline = "Submitter-processed"
    result.save()

    # Create a ComputedFile for the original file,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = job_context["output_file_path"]
        computed_file.filename = job_context['output_file_path'].split('/')[-1]
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        # computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
        # TODO here: delete local file after S3 sync
        computed_file.save()
    except Exception:
        logger.error("Exception caught while moving file %s",
                         raw_path,
                         processor_job=job_context["job_id"])

        failure_reason = "Exception caught while moving file {}".format(file.name)
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

def no_op_processor(job_id: int) -> None:
    return utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _convert_genes,
                        _create_result,
                        utils.end_job])

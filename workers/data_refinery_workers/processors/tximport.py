import boto3
import glob
import io
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import tarfile

from botocore.client import Config
from django.conf import settings
from django.db import transaction
from django.utils import timezone
from typing import Dict, List
import numpy as np
import pandas as pd

from data_refinery_common.job_lookup import Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Experiment,
    ExperimentSampleAssociation,
    OrganismIndex,
    Pipeline,
    Processor,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils, salmon

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client('s3', config=Config(signature_version='s3v4'))
logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "processor_job_"
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

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

    # Create a directory specific to this processor job.
    # (A single sample could belong to multiple experiments, meaning
    # that it could be run more than once, potentially even at the
    # same time.)
    job_context["work_dir"] = os.path.join(LOCAL_ROOT_DIR,
                                           job_context["job_dir_prefix"]) + "/"
    job_context["temp_dir"] = job_context["work_dir"] + "temp/"
    os.makedirs(job_context["work_dir"], exist_ok=True)
    os.makedirs(job_context["temp_dir"], exist_ok=True)

    # Technically unsafe, but if either of these objects don't exist we need to fail anyway.
    sample = job_context["job"].original_files.first().samples.first()
    job_context['sample'] = sample
    job_context['samples'] = []
    job_context['organism'] = sample.organism
    job_context["success"] = True
    job_context["is_tximport_only"] = True

    job_context["computed_files"] = []
    job_context["smashable_files"] = []

    archive_file = ComputedFile.objects.filter(
        result_id__in=sample.results.values('id'),
        filename__startswith='result-',
        filename__endswith='.tar.gz'
    ).latest('last_modified')

    archive_path = archive_file.get_synced_file_path()

    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(job_context["temp_dir"])

    return job_context


def tximport(job_id: int) -> None:
    """Main processor function for the Tximport Processor.

    Runs salmon quant command line tool, specifying either a long or
    short read length. Also runs FastQC, MultiQC, and Salmontools.
    """
    pipeline = Pipeline(name=utils.PipelineEnum.TXIMPORT.value)
    final_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                       [utils.start_job,
                        _set_job_prefix,
                        _prepare_files,
                        salmon._find_or_download_index,
                        salmon._tximport,
                        utils.end_job])
    return final_context

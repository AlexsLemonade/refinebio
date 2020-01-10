import os
from typing import Dict

from django.conf import settings
from django.db import transaction
from django.utils import timezone

import boto3
import numpy as np
import pandas as pd
from botocore.client import Config

from data_refinery_common.job_lookup import Downloaders, PipelineEnum
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
from data_refinery_workers.processors import salmon, utils

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))
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
    """
    logger.debug("Preparing files..")

    # Create a directory specific to this processor job.
    # (A single sample could belong to multiple experiments, meaning
    # that it could be run more than once, potentially even at the
    # same time.)
    job_context["work_dir"] = os.path.join(LOCAL_ROOT_DIR, job_context["job_dir_prefix"]) + "/"
    os.makedirs(job_context["work_dir"], exist_ok=True)

    # Technically unsafe, but if either of these objects don't exist we need to fail anyway.
    sample = job_context["job"].original_files.first().samples.first()
    job_context["sample"] = sample
    job_context["samples"] = []
    job_context["organism"] = sample.organism
    job_context["success"] = True
    job_context["is_tximport_only"] = True

    job_context["computed_files"] = []
    job_context["smashable_files"] = []

    return job_context


def tximport(job_id: int) -> None:
    """Main processor function for the Tximport Processor.

    Runs tximport command line tool on an experiment.
    """
    pipeline = Pipeline(name=PipelineEnum.TXIMPORT.value)
    final_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _set_job_prefix,
            _prepare_files,
            salmon.get_tximport_inputs,
            salmon._find_or_download_index,
            salmon.tximport,
            utils.end_job,
        ],
    )
    return final_context

"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals

from enum import Enum

from django.conf import settings

import boto3
import nomad
from nomad.api.exceptions import URLNotFoundNomadException

from data_refinery_common.job_lookup import (
    SMASHER_JOB_TYPES,
    Downloaders,
    ProcessorPipeline,
    SurveyJobTypes,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import DownloaderJob, ProcessorJob, SurveyJob
from data_refinery_common.utils import get_env_variable, get_volume_index

logger = get_and_configure_logger(__name__)


AWS_REGION = get_env_variable(
    "AWS_REGION", "us-east-1"
)  # Default to us-east-1 if the region variable can't be found
AWS_BATCH_QUEUE_NAME = get_env_variable("REFINEBIO_JOB_QUEUE_NAME")

# These two constants refer to image names that should be used for
# multiple jobs.
NOMAD_TRANSCRIPTOME_JOB = "TRANSCRIPTOME_INDEX"
NOMAD_DOWNLOADER_JOB = "DOWNLOADER"
NONE_JOB_ERROR_TEMPLATE = "send_job was called with NONE job_type: {} for {} job {}"


def get_job_name(job_type, job_id):
    if (
        job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG
        or job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT
    ):
        return NOMAD_TRANSCRIPTOME_JOB
    elif job_type is ProcessorPipeline.SALMON or job_type is ProcessorPipeline.TXIMPORT:
        # Tximport uses the same job specification as Salmon.
        return ProcessorPipeline.SALMON.value
    elif job_type is ProcessorPipeline.AFFY_TO_PCL:
        return ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type is ProcessorPipeline.NO_OP:
        return ProcessorPipeline.NO_OP.value
    elif job_type is ProcessorPipeline.ILLUMINA_TO_PCL:
        return ProcessorPipeline.ILLUMINA_TO_PCL.value
    elif job_type is ProcessorPipeline.SMASHER:
        return ProcessorPipeline.SMASHER.value
    elif job_type is ProcessorPipeline.JANITOR:
        return ProcessorPipeline.JANITOR.value
    elif job_type is ProcessorPipeline.QN_REFERENCE:
        return ProcessorPipeline.QN_REFERENCE.value
    elif job_type is ProcessorPipeline.CREATE_COMPENDIA:
        return ProcessorPipeline.CREATE_COMPENDIA.value
    elif job_type is ProcessorPipeline.CREATE_QUANTPENDIA:
        return ProcessorPipeline.CREATE_QUANTPENDIA.value
    elif job_type is ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL:
        # Agilent twocolor uses the same job specification as Affy.
        return ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type in list(Downloaders):
        return NOMAD_DOWNLOADER_JOB
    elif job_type in list(SurveyJobTypes):
        return job_type.value
    elif job_type is Downloaders.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Downloader", job_id))
    elif job_type is ProcessorPipeline.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Processor", job_id))
    else:
        raise ValueError("Invalid job_type: {}".format(job_type.value))


def is_job_processor(job_type):
    if job_type in list(Downloaders):
        return False
    elif job_type in list(SurveyJobTypes):
        return False

    return True


def send_job(job_type: Enum, job, is_dispatch=False) -> bool:
    job_name = get_job_name(job_type, job.id)
    is_processor = is_job_processor(job_type)

    if settings.AUTO_DISPATCH_NOMAD_JOBS:
        # We only want to dispatch processor jobs directly.
        # Everything else will be handled by the Foreman, which will increment the retry counter.
        should_dispatch = is_processor or is_dispatch or (not settings.RUNNING_IN_CLOUD)
    else:
        should_dispatch = is_dispatch  # only dispatch when specifically requested to

    # Temporary until the foreman will dispatch these correctly.
    should_dispatch = True

    if should_dispatch:
        batch = boto3.client("batch", region_name=AWS_REGION)

        # Smasher doesn't need to be on a specific instance since it will
        # download all the data to its instance anyway.
        if isinstance(job, ProcessorJob) and job_type not in SMASHER_JOB_TYPES:
            # Make sure this job goes to the correct EBS resource.
            # If this is being dispatched for the first time, make sure that
            # we store the currently attached index.
            # If this is being dispatched by the Foreman, it should already
            # have an attached volume index, so use that.
            if job.volume_index is None:
                job.volume_index = get_volume_index()
                job.save()
            job_name = job_name + "_" + job.volume_index + "_" + str(job.ram_amount)
        elif isinstance(job, SurveyJob):
            job_name = job_name + "_" + str(job.ram_amount)
        elif isinstance(job, DownloaderJob):
            # volume_index = job.volume_index if settings.RUNNING_IN_CLOUD else "0"
            volume_index = job.volume_index or "0"
            job_name = job_name + "_" + volume_index + "_" + str(job.ram_amount)

        try:
            batch_response = batch.submit_job(
                jobName=job_name + f"_{job.id}",
                jobQueue=AWS_BATCH_QUEUE_NAME,
                jobDefinition=job_name,
                parameters={"job_name": job_type.value, "job_id": str(job.id)},
            )
            job.nomad_job_id = batch_response["jobId"]
            job.save()
            return True
        except Exception as e:
            logger.info(
                "Unable to Dispatch Batch Job.",
                job_name=job_type.value,
                job_id=str(job.id),
                reason=str(e),
            )
            raise
    else:
        job.num_retries = job.num_retries - 1
        job.save()
    return True

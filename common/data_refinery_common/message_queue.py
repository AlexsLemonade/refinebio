"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals

from enum import Enum

from django.conf import settings

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
from data_refinery_common.utils import (
    get_env_variable,
    get_env_variable_gracefully,
    get_volume_index,
)

logger = get_and_configure_logger(__name__)


# These two constants refer to image names that should be used for
# multiple jobs.
NOMAD_TRANSCRIPTOME_JOB = "TRANSCRIPTOME_INDEX"
NOMAD_DOWNLOADER_JOB = "DOWNLOADER"
NONE_JOB_ERROR_TEMPLATE = "send_job was called with NONE job_type: {} for {} job {}"


def send_job(job_type: Enum, job, is_dispatch=False) -> bool:
    """Queues a worker job by sending a Nomad Job dispatch message.

    job_type must be a valid Enum for ProcessorPipelines or
    Downloaders as defined in data_refinery_common.job_lookup.
    job must be an existing ProcessorJob or DownloaderJob record.

    Returns True if the job was successfully dispatch, return False otherwise.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

    is_processor = True
    if (
        job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG
        or job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT
    ):
        nomad_job = NOMAD_TRANSCRIPTOME_JOB
    elif job_type is ProcessorPipeline.SALMON or job_type is ProcessorPipeline.TXIMPORT:
        # Tximport uses the same job specification as Salmon.
        nomad_job = ProcessorPipeline.SALMON.value
    elif job_type is ProcessorPipeline.AFFY_TO_PCL:
        nomad_job = ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type is ProcessorPipeline.NO_OP:
        nomad_job = ProcessorPipeline.NO_OP.value
    elif job_type is ProcessorPipeline.ILLUMINA_TO_PCL:
        nomad_job = ProcessorPipeline.ILLUMINA_TO_PCL.value
    elif job_type is ProcessorPipeline.SMASHER:
        nomad_job = ProcessorPipeline.SMASHER.value
    elif job_type is ProcessorPipeline.JANITOR:
        nomad_job = ProcessorPipeline.JANITOR.value
    elif job_type is ProcessorPipeline.QN_REFERENCE:
        nomad_job = ProcessorPipeline.QN_REFERENCE.value
    elif job_type is ProcessorPipeline.CREATE_COMPENDIA:
        nomad_job = ProcessorPipeline.CREATE_COMPENDIA.value
    elif job_type is ProcessorPipeline.CREATE_QUANTPENDIA:
        nomad_job = ProcessorPipeline.CREATE_QUANTPENDIA.value
    elif job_type is ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL:
        # Agilent twocolor uses the same job specification as Affy.
        nomad_job = ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type in list(Downloaders):
        nomad_job = NOMAD_DOWNLOADER_JOB
        is_processor = False
    elif job_type in list(SurveyJobTypes):
        nomad_job = job_type.value
        is_processor = False
    elif job_type is Downloaders.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Downloader", job_id))
    elif job_type is ProcessorPipeline.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Processor", job_id))
    else:
        raise ValueError("Invalid job_type: {}".format(job_type.value))

    logger.debug(
        "Queuing %s nomad job to run job %s with id %d.", nomad_job, job_type.value, job.id
    )

    # We only want to dispatch processor jobs directly.
    # Everything else will be handled by the Foreman, which will increment the retry counter.
    if is_processor or is_dispatch or (not settings.RUNNING_IN_CLOUD):

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
            nomad_job = nomad_job + "_" + job.volume_index + "_" + str(job.ram_amount)
        elif isinstance(job, SurveyJob):
            nomad_job = nomad_job + "_" + str(job.ram_amount)
        elif isinstance(job, DownloaderJob):
            volume_index = job.volume_index if settings.RUNNING_IN_CLOUD else "0"
            nomad_job = nomad_job + "_" + volume_index + "_" + str(job.ram_amount)

        try:
            nomad_response = nomad_client.job.dispatch_job(
                nomad_job, meta={"JOB_NAME": job_type.value, "JOB_ID": str(job.id)}
            )
            job.nomad_job_id = nomad_response["DispatchedJobID"]
            job.save()
            return True
        except URLNotFoundNomadException:
            logger.info(
                "Dispatching Nomad job of type %s for job spec %s to host %s and port %s failed.",
                job_type,
                nomad_job,
                nomad_host,
                nomad_port,
                job=str(job.id),
            )
            raise
        except Exception as e:
            logger.info(
                "Unable to Dispatch Nomad Job.",
                job_name=job_type.value,
                job_id=str(job.id),
                reason=str(e),
            )
            raise
    else:
        job.num_retries = job.num_retries - 1
        job.save()
    return True

"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals
from enum import Enum
import nomad
from nomad.api.exceptions import URLNotFoundNomadException
from data_refinery_common.utils import get_env_variable
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


# These two constants refer to image names that should be used for
# multiple jobs.
NOMAD_TRANSCRIPTOME_JOB = "TRANSCRIPTOME_INDEX"
NOMAD_DOWNLOADER_JOB = "DOWNLOADER"
NONE_JOB_ERROR_TEMPLATE = "send_job was called with NONE job_type: {} for {} job {}"


def send_job(job_type: Enum, job) -> None:
    """Queues a worker job by sending a Nomad Job dispatch message.

    job_type must be a valid Enum for ProcessorPipelines or
    Downloaders as defined in data_refinery_common.job_lookup.
    job must be an existing ProcessorJob or DownloaderJob record.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=5)

    if job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG \
       or job_type is ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT:
        nomad_job = NOMAD_TRANSCRIPTOME_JOB
    elif job_type is ProcessorPipeline.SALMON:
        nomad_job = ProcessorPipeline.SALMON.value
    elif job_type is ProcessorPipeline.AFFY_TO_PCL:
        nomad_job = ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type is ProcessorPipeline.NO_OP:
        nomad_job = ProcessorPipeline.NO_OP.value
    elif job_type is ProcessorPipeline.ILLUMINA_TO_PCL:
        nomad_job = ProcessorPipeline.ILLUMINA_TO_PCL.value
    elif job_type is ProcessorPipeline.SMASHER:
        nomad_job = ProcessorPipeline.SMASHER.value
    elif job_type is ProcessorPipeline.AGILENT_TWOCOLOR_TO_PCL:
        # Agilent twocolor uses the same job specification as Affy.
        nomad_job = ProcessorPipeline.AFFY_TO_PCL.value
    elif job_type is Downloaders.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Downloader", job_id))
    elif job_type is ProcessorPipeline.NONE:
        logger.warn("Not queuing %s job.", job_type, job_id=job_id)
        raise ValueError(NONE_JOB_ERROR_TEMPLATE.format(job_type.value, "Processor", job_id))
    elif job_type in list(Downloaders):
        nomad_job = NOMAD_DOWNLOADER_JOB
    else:
        raise ValueError("Invalid job_type: {}".format(job_type.value))

    logger.info("Queuing %s nomad job to run job %s with id %d.",
                nomad_job,
                job_type.value,
                job.id)
    try:
        nomad_response = nomad_client.job.dispatch_job(nomad_job, meta={"JOB_NAME": job_type.value,
                                                                        "JOB_ID": str(job.id)})
        job.nomad_job_id = nomad_response["DispatchedJobID"]
        job.save()
    except URLNotFoundNomadException:
        logger.error("Dispatching Nomad job of type %s to host %s and port %s failed.",
                     job_type, nomad_host, nomad_port, job=str(job))
        raise

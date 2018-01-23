"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals
from enum import Enum
import nomad
from data_refinery_common.utils import get_env_variable
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


# There are currently two Nomad Job Specifications defined in
# workers/downloader.nomad.tpl and workers/processor.nomad.tpl.
# These constants are the identifiers for those two job specifications.
NOMAD_PROCESSOR_JOB = "PROCESSOR"
NOMAD_DOWNLOADER_JOB = "DOWNLOADER"


def send_job(job_type: Enum, job_id: int) -> None:
    """Queues a worker job by sending a Nomad Job dispatch message.

    job_type must be a valid Enum for ProcessorPipelines or
    Downloaders as defined in data_refinery_common.job_lookup.
    job_id must correspond to an existing ProcessorJob or
    DownloaderJob record.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_client = nomad.Nomad(nomad_host, timeout=5)

    # Once I have every job specced out with its own Nomad job, this
    # code can change and the meta won't need "JOB_NAME" in it because
    # the just specifying the nomad_job to dispatch will be enough.
    if job_type in list(ProcessorPipeline):
        nomad_job = NOMAD_PROCESSOR_JOB
    elif job_type in list(Downloaders):
        nomad_job = NOMAD_DOWNLOADER_JOB
    else:
        raise ValueError("Invalid job_type.")

    logger.info("Queuing %s nomad job to run DR job %s with id %d.",
                nomad_job,
                job_type.value,
                job_id)
    nomad_client.job.dispatch_job(nomad_job, meta={"JOB_NAME": job_type.value,
                                                   "JOB_ID": str(job_id)})

    # This is what that will look like eventually:
    # nomad_client.job.dispatch_job(job_type.value, meta={"JOB_ID": str(job_id)})

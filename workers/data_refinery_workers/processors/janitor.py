import os
import random
import shutil
import string
import subprocess
import time
import warnings

from django.utils import timezone
from nomad import Nomad
from nomad.api.exceptions import BaseNomadException, URLNotFoundNomadException
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Pipeline,
    Processor,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    ProcessorJob
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
logger = get_and_configure_logger(__name__)


def _find_and_remove_expired_jobs(job_context):
    """ Finds expired jobs and removes their working directories """

    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=15)

    job_context['deleted_items'] = []

    for item in os.listdir(LOCAL_ROOT_DIR):

        # Processor job working directories
        if 'processor_job_' in item:

            # TX Index jobs are the only ones who are allowed to hang around
            # after their jobs are finished. They're marked with an _index in their path.
            if '_index' in item:
                continue

            job_id = item.split('processor_job_')[1]

            # Okay, does this job exist?
            try:
                job = ProcessorJob.objects.get(id=job_id)

                # Is this job running?
                try:
                    job_status = nomad_client.job.get_job(job.nomad_job_id)["Status"]

                    # This job is running, don't delete  the working directory.
                    if job_status == "running":
                        continue
                    if job_status == "pending":
                        continue
                except URLNotFoundNomadException as e:
                    # Nomad has no record of this job, meaning it has likely been GC'd after death.
                    # It can be purged.
                    pass
                except BaseNomadException as e:
                    # If we can't currently access Nomad,
                    # just continue until we can again.
                    continue
                except Exception as e:
                    # This job is likely vanished.
                    # Or, possibly, another Nomad error outside of BaseNomadException.
                    continue
            except ProcessorJob.DoesNotExist:
                # This job has vanished from the DB - clean it up!
                logger.error("Janitor found no record of " + item + " - why?")
                pass
            except Exception:
                # We're unable to connect to the DB right now (or something), so hold onto it for right now.
                logger.exception("Problem finding job record for " + item + " - why?")
                continue

            # Delete it!
            try:
                to_delete = LOCAL_ROOT_DIR + '/' + item
                logger.info("Janitor deleting " + to_delete, contents=str(os.listdir(to_delete)), job_id=job_id)
                shutil.rmtree(to_delete)
                job_context['deleted_items'].append(to_delete)
            except Exception as e:
                # This job is likely vanished. No need for this directory.
                pass

        # There may be successful processors
        if 'SRP' in item or 'ERP' in item or 'DRR' in item:
            sub_path = os.path.join(LOCAL_ROOT_DIR, item)
            for sub_item in os.listdir(sub_path):
                try:
                    sample = Sample.objects.get(accession_code=sub_item)
                    if sample.computed_files.count() == 0:
                        # This doesn't have any associated computed files - leave it be.
                        continue
                except Sample.DoesNotExist:
                    # Interesting. This shouldn't happen at all.
                    continue
                except Exception:
                    # We can't contact the DB right now, skip deletion.
                    continue

                try:
                    sub_item_path = os.path.join(sub_path, sub_item)
                    logger.info("Janitor deleting " + sub_item_path, contents=str(os.listdir(sub_item_path)))
                    shutil.rmtree(sub_item_path)
                    job_context['deleted_items'].append(sub_item_path)
                except Exception as e:
                    # This job is likely vanished. No need for this directory.
                    pass

    job_context['success'] = True
    return job_context

def run_janitor(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.JANITOR.value)
    job_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                       [utils.start_job,
                        _find_and_remove_expired_jobs,
                        utils.end_job])
    return job_context

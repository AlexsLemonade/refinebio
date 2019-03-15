"""
This command will run the Foreman's main function monitor_jobs.
This will cause the Foreman to check for a number of different
failures for both the DownloaderJobs and ProcessorJobs and requeue
those jobs it detects as failed.
"""

import random
import sys
import time
from typing import Dict, List

from django.core.management.base import BaseCommand
from nomad import Nomad

from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
)
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.utils import get_env_variable, get_active_volumes


logger = get_and_configure_logger(__name__)


MAX_JOBS_FOR_THIS_MODE = 10000


def build_completion_list(organism: Organism) -> List[Dict]:
    """Builds a list of stats on how close to completion experiments are.

    For all experiments related to organism, adds a dictionary to the
    list which contains the experiment, the number of processed
    samples, and the number of unprocessed samples under the keys:
    experiment, processed, and unprocessed respectively.
    """
    experiments = organism.experiments.filter(technology='RNA-SEQ')

    completion_list = []
    for experiment in experiments:
        processed_samples = set()
        unprocessed_samples = set()

        for sample in experiment.samples.all():
            if sample.is_processed:
                continue

            sample_is_processed = False
            for original_file in sample.original_files.all():
                if sample_is_processed:
                    break

                for processor_job in original_file.processor_jobs.all():
                    if processor_job.success:
                        sample_is_processed = True
                        break

            if sample_is_processed:
                processed_samples.add(sample)
            else:
                unprocessed_samples.add(sample)

        completion_list.append({
            "experiment": experiment,
            "unprocessed": unprocessed_samples,
            "processed": processed_samples
        })

    return completion_list


def build_prioritized_jobs_list(organism: Organism) -> List:
    """Returns a prioritized list of stats on how close to completion experiments are.

    The list is sorted from highest completion percent to lowest and
    has the same structure as the list returned by
    build_completion_list.
    """
    completion_list = build_completion_list(organism)

    def calculate_completion_percentage(experiment_stats_dict):
        num_processed = len(experiment_stats_dict["processed"])
        num_unprocessed = len(experiment_stats_dict["unprocessed"])

        if num_processed + num_unprocessed == 0:
            return 0

        return num_processed / (num_processed + num_unprocessed)

    completion_list.sort(reverse=True, key=calculate_completion_percentage)

    prioritized_job_list = []
    for experiment_stats_dict in completion_list:
        unprocessed_samples = experiment_stats_dict["unprocessed"]
        logger.info(
            "Experiment %s has %d / %d samples processed.",
            experiment_stats_dict["experiment"],
            len(unprocessed_samples),
            (len(unprocessed_samples) + len(experiment_stats_dict["processed"]))
        )

        for sample in unprocessed_samples:
            processor_jobs = list(sample.get_processor_jobs())
            downloader_jobs = list(sample.get_downloader_jobs())
            if processor_jobs:
                # We want to requeue the most recently created
                # processor job, so sort by id since they are
                # autoincrementing.
                processor_jobs.sort(
                    reverse=True,
                    key=lambda x: x.id
                )
                prioritized_job_list.append(processor_jobs[0])
            elif downloader_jobs:
                downloader_jobs.sort(
                    reverse=True,
                    key=lambda x: x.id
                )
                prioritized_job_list.append(downloader_jobs[0])

    return prioritized_job_list


def requeue_job(job):
    """Requeues a job regardless of whether it is a DownloaderJob or ProcessorJob.

    This function reuses a lot of logic from requeue_downloader_job
    and requeue_processor_job from the main namespace for the
    Foreman. However there's additional logic there that we don't
    want, because we explicitly want to requeue these jobs regardless
    of how many times they've been retried.
    """
    # All new jobs are going to be set at 2 retries so they only get
    # tried once. Presumably all of these have failed at least once
    # already, so there may be a good reason.
    num_retries = 1
    if isinstance(job, ProcessorJob):
        # We don't want these jobs to sit in our queue because the
        # volume we assigned isn't available, so only use active
        # volumes. Also in order to spread them around do so randomly.
        volume_index = random.choice(list(get_active_volumes()))
        new_job = ProcessorJob(
            num_retries=num_retries,
            pipeline_applied=job.pipeline_applied,
            ram_amount=job.ram_amount,
            volume_index=volume_index
        )
        new_job.save()

        for original_file in job.original_files.all():
            ProcessorJobOriginalFileAssociation.objects.get_or_create(
                processor_job=new_job,
                original_file=original_file
            )

        job_type = ProcessorPipeline[job.pipeline_applied]
    elif isinstance(job, DownloaderJob):
        new_job = DownloaderJob(
            num_retries=num_retries,
            downloader_task=job.downloader_task,
            accession_code=job.accession_code
        )
        new_job.save()

        for original_file in job.original_files.all():
            DownloaderJobOriginalFileAssociation.objects.get_or_create(
                downloader_job=new_job,
                original_file=original_file
            )

        job_type = Downloaders[job.downloader_task]
    else:
        raise ValueError("Told to requeue a job that's not a ProcessorJob nor DownloaderJob!")

    try:
        if send_job(job_type, job=new_job, is_dispatch=True):
            job.retried = True
            job.success = False
            job.retried_job = new_job
            job.save()
        else:
            logger.error(
                ("Failed to requeue %s which had ID %d with a new %s "
                 "with ID %d because send_job returned false."),
                type(job).__name__,
                job.id,
                type(job).__name__,
                new_job.id
            )
            # Can't communicate with nomad just now, leave the job for a later loop.
            new_job.delete()
            return False
    except:
        logger.exception(
            ("Failed to requeue %s which had ID %d with a new %s "
             "with ID %d because send_job raised an exception."),
            type(job).__name__,
            job.id,
            type(job).__name__,
            new_job.id
        )
        # Can't communicate with nomad just now, leave the job for a later loop.
        new_job.delete()
        return False

    return True


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--organism-name",
            type=str,
            help=("The name of the organism you would like to shepherd to completion. "
                  "Must be in all caps and any spaces should be replaced with an underscore. "
                  "I should make this less strict.")
        )

    def handle(self, *args, **options):
        """Requeues all unprocessed RNA-Seq samples for an organism.
        """
        if options["organism_name"] is None:
            logger.error("You must specify an organism-name.")
            sys.exit(1)
        else:
            organism_name = options["organism_name"]

        organism = Organism.objects.get(name=organism_name)

        prioritized_job_list = build_prioritized_jobs_list(organism)

        if not len(prioritized_job_list):
            logger.info("Found no samples that need to be processed. I guess I'm done!")
            sys.exit(0)

        logger.info(
            "Found %d samples that need to be processed. Beginning to queue jobs!",
            len(prioritized_job_list)
        )

        nomad_host = get_env_variable("NOMAD_HOST")
        nomad_port = get_env_variable("NOMAD_PORT", "4646")
        nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

        while(len(prioritized_job_list) > 0):
            len_all_jobs = len(nomad_client.jobs.get_jobs())

            num_short_from_max = MAX_JOBS_FOR_THIS_MODE - len_all_jobs
            if num_short_from_max > 0:
                for i in range(num_short_from_max):
                    if len(prioritized_job_list) > 0:
                        requeue_job(prioritized_job_list.pop(0))

            # Wait 10 minutes in between queuing additional work to
            # give it time to actually get done.
            if len(prioritized_job_list) > 0:
                logger.info("Sleeping for 5 minutes while jobs get done.")
                time.sleep(600)

        logger.info("Successfully requeued all jobs for unprocessed %s samples.", organism_name)

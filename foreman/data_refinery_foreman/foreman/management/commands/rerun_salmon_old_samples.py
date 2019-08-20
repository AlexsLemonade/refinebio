""" This management command will search for experiments where samples haven't been processed and
that multiple versions of salmon have been used on the samples.
"""

import random
import sys
import time
from typing import Dict, List

from django.core.management.base import BaseCommand
from nomad import Nomad
from django.db.models import OuterRef, Subquery

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
    ComputationalResult,
)
from data_refinery_common.job_lookup import ProcessorEnum, ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.rna_seq import get_quant_results_for_experiment
from data_refinery_common.utils import get_env_variable, get_active_volumes
from data_refinery_common.job_management import create_processor_job_for_original_files

logger = get_and_configure_logger(__name__)

def update_salmon_versions(experiment: Experiment):
    quant_results = get_quant_results_for_experiment(experiment)
    salmon_versions = list(quant_results.order_by('-organism_index__created_at')\
                                   .values_list('organism_index__salmon_version', flat=True)\
                                   .distinct())
    
    if len(salmon_versions) <= 1:
        # only apply this command on experiments that have more than one salmon version applied on their samples
        return

    latest_salmon_version = salmon_versions[0]

    # find the samples that were not processed with `latest_salmon_version` and trigger new processor jobs for them
    newest_computational_results = ComputationalResult.objects.all()\
        .filter(
            samples=OuterRef('id'),
            processor__name=ProcessorEnum.SALMON_QUANT.value['name']
        )\
        .order_by('-created_at')
    
    samples = experiment.samples.all().annotate(
            salmon_version=Subquery(newest_computational_results.values('organism_index__salmon_version')[:1])
        )\
        .exclude(salmon_version=latest_salmon_version)

    # create new processor jobs for the samples that were run with an older salmon version
    for sample in samples:
        original_files = list(sample.original_files.all())
        create_processor_job_for_original_files(original_files)

def update_salmon_all_experiments():
    """Creates a tximport job for all eligible experiments."""
    eligible_experiments = Experiment.objects.all()\
        .filter(technology='RNA-SEQ', num_processed_samples=0)\
    
    for experiment in eligible_experiments:
        update_salmon_versions(experiment)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession-code",
            type=str,
            help=("Optional parameter to allow running the command on a single experiment")
        )

    def handle(self, *args, **options):
        if options["accession_code"] is None:
            update_salmon_all_experiments()
        else:
            experiment = Experiment.objects.all().filter(accession_code=options["accession_code"]).first()
            if not experiment:
                logger.error("The provided experiment accession code was not valid.")
                sys.exit(1)
            else:
                update_salmon_versions(experiment)

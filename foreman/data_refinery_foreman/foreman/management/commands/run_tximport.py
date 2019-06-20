"""This management command will queue tximport jobs for experiments
which haven't been completely processed yet. It will only do so for
experiments which are over the minimum thresholds for running tximport
and which haven't had tximport run yet. It also will only run for
experiments with a single organism since we don't yet have logic to
split experiments into multiple tximport jobs.
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
from data_refinery_common.rna_seq import get_quant_results_for_experiment, should_run_tximport
from data_refinery_common.utils import get_env_variable, get_active_volumes


logger = get_and_configure_logger(__name__)

class Command(BaseCommand):

    def handle(self, *args, **options):

        eligible_experiments = Experiment.objects.extra(
            select={'num_organisms':'cardinality(organism_names)'}
        ).filter(
            'num_organisms'=1,
            technology='RNA-SEQ',
            num_processed_samples=0
        ).prefetch_related('samples__results__computed_files')

        # Next is to figure out how many samples were processed for
        # each experiment. Should be able to reuse code from salmon
        # cause it does this stuff.
        tximport_pipeline = ProcessorPipeline.TXIMPORT
        for experiment in eligible_experiments:
            quant_results = get_quant_results_for_experiment(experiment)
            if should_run_tximport(experiment, len(quant_results), True):
                processor_job = ProcessorJob()
                processor_job.pipeline_applied = tximport_pipeline.value
                processor_job.ram_amount = 8192
                processor_job.save()

                assoc = ProcessorJobOriginalFileAssociation()
                # Any original file linked to any sample of the
                # experiment will work. Tximport is somewhat special
                # in that it doesn't actuallhy use original files so
                # this is just used to point to the experiment.
                assoc.original_file = experiment.samples.all()[0].original_files.all()[0]
                assoc.processor_job = processor_job
                assoc.save()

                try:
                    send_job(tximport_pipeline, processor_job)
                except:
                    # If we cannot queue the job now the Foreman will do
                    # it later.
                    pass

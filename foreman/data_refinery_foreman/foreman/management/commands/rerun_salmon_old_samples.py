""" This management command will search for experiments where samples haven't been processed and
that multiple versions of salmon have been used on the samples.
"""

import random
import sys

from django.core.management.base import BaseCommand
from nomad import Nomad
from django.db.models import Count
from django.db.models.expressions import Q

from data_refinery_common.models import (
    Experiment,
    ProcessorJob,
)
from data_refinery_common.job_lookup import ProcessorEnum, ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.rna_seq import get_quant_results_for_experiment
from data_refinery_common.job_management import create_downloader_job
from data_refinery_foreman.foreman import main

logger = get_and_configure_logger(__name__)


def update_salmon_versions(experiment: Experiment):
    quant_results = (
        get_quant_results_for_experiment(experiment, filter_old_versions=False)
        .order_by("-organism_index__created_at")
        .prefetch_related("organism_index")
        .prefetch_related("samples__original_files")
    )

    total_samples_queued = 0
    latest_salmon_version = None
    for quant_result in quant_results:
        if not latest_salmon_version:
            # we can safely ignore the latest salmon version, that will be the first
            # quant result. Note we are ordering by -organism_index__created_at
            latest_salmon_version = quant_result.organism_index.salmon_version
        elif latest_salmon_version != quant_result.organism_index.salmon_version:
            # we found a quant result associated with an experiment where we need to run salmon
            # hopefully each computational result is associated with a single sample
            for sample in quant_result.samples.all():
                original_files = list(sample.original_files.all())

                if not len(original_files):
                    continue

                # Ensure that there's no processor jobs for these original files that the foreman
                # might want to retry (failed | hung | lost)
                has_open_processor_job = (
                    ProcessorJob.objects.filter(
                        original_files=original_files[0],
                        pipeline_applied=ProcessorPipeline.SALMON,
                    )
                    .filter(
                        Q(success=False, retried=False, no_retry=False)
                        | Q(
                            success=None,
                            retried=False,
                            no_retry=False,
                            start_time__isnull=False,
                            end_time=None,
                            nomad_job_id__isnull=False,
                        )
                        | Q(
                            success=None,
                            retried=False,
                            no_retry=False,
                            start_time=None,
                            end_time=None,
                        )
                    )
                    .exists()
                )
                if has_open_processor_job:
                    continue

                create_downloader_job(original_files, force=True)
                total_samples_queued += 1

    logger.info(
        "Re-ran Salmon for %d samples in experiment %s.",
        total_samples_queued,
        experiment.accession_code,
    )


def update_salmon_all_experiments():
    """Creates a tximport job for all eligible experiments."""
    eligible_experiments = (
        Experiment.objects.filter(technology="RNA-SEQ", num_processed_samples=0)
        .annotate(
            num_salmon_versions=Count(
                "samples__results__organism_index__salmon_version",
                distinct=True,
                filter=Q(
                    samples__results__processor__name=ProcessorEnum.SALMON_QUANT.value[
                        "name"
                    ]
                ),
            )
        )
        .filter(num_salmon_versions__gt=1)
    )

    for experiment in eligible_experiments:
        update_salmon_versions(experiment)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession-code",
            type=str,
            help=(
                "Optional parameter to allow running the command on a single experiment"
            ),
        )

    def handle(self, *args, **options):
        if options["accession_code"] is None:
            update_salmon_all_experiments()
        else:
            experiment = Experiment.objects.filter(
                accession_code=options["accession_code"]
            ).first()
            if not experiment:
                logger.error("The provided experiment accession code was not valid.")
                sys.exit(1)
            else:
                update_salmon_versions(experiment)

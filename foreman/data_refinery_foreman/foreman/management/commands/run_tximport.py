"""This management command will queue tximport jobs for experiments
which haven't been completely processed yet. It will only do so for
experiments which are over the minimum thresholds for running tximport
and which haven't had tximport run yet. It also will only run for
experiments with a single organism since we don't yet have logic to
split experiments into multiple tximport jobs.
"""

from django.core.management.base import BaseCommand
from django.db.models import Count

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    Experiment,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
)
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_common.rna_seq import get_quant_results_for_experiment, should_run_tximport

logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000


def run_tximport():
    """Creates a tximport job for all eligible experiments."""
    eligible_experiments = (
        Experiment.objects.annotate(num_organisms=Count("organisms"))
        .filter(num_organisms=1, technology="RNA-SEQ", num_processed_samples=0)
        .prefetch_related("samples__results")
    )

    paginator = Paginator(eligible_experiments, PAGE_SIZE)
    page = paginator.page()

    # Next is to figure out how many samples were processed for
    # each experiment. Should be able to reuse code from salmon
    # cause it does this stuff.
    tximport_pipeline = ProcessorPipeline.TXIMPORT

    created_jobs = []

    while True:
        creation_count = 0

        for experiment in page.object_list:
            quant_results = get_quant_results_for_experiment(experiment)

            if should_run_tximport(experiment, quant_results, True):
                processor_job = ProcessorJob()
                processor_job.pipeline_applied = tximport_pipeline.value
                processor_job.ram_amount = 32768
                processor_job.save()

                assoc = ProcessorJobOriginalFileAssociation()
                # Any original file linked to any sample of the
                # experiment will work. Tximport is somewhat special
                # in that it doesn't actuallhy use original files so
                # this is just used to point to the experiment.
                assoc.original_file = experiment.samples.all()[0].original_files.all()[0]
                assoc.processor_job = processor_job
                assoc.save()

                creation_count += 1
                created_jobs.append(processor_job)

                try:
                    send_job(tximport_pipeline, processor_job)
                except Exception:
                    # If we cannot queue the job now the Foreman will do
                    # it later.
                    pass

        logger.info("Created %d tximport jobs for experiments past the thresholds.", creation_count)

        if not page.has_next():
            break
        else:
            page = paginator.page(page.next_page_number())

    return created_jobs


class Command(BaseCommand):
    def handle(self, *args, **options):
        """This command just calls run_tximport.

        The functionality has been broken out into a separate function
        to make testing easy.
        """
        run_tximport()

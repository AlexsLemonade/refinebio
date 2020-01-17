import sys

from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    Experiment,
    OriginalFile,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
)
from data_refinery_workers.processors import tximport

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--experiment-accession",
            type=str,
            help=("Accession code for the experiment to run tximport on."),
        )

    def handle(self, *args, **options):
        """Run tximport on an experiment for the samples that have been processed by Salmon.
        """
        if options["experiment_accession"] is None:
            logger.error("The --experiment-accession argument must be provided")
            sys.exit(1)
        else:
            accession_code = options["experiment_accession"]

        # Find an OriginalFile associated with one of the Experiment's
        # samples which had Salmon run successfully on it.
        original_file = None
        experiment = Experiment.objects.get(accession_code=accession_code)
        for sample in experiment.samples.all():
            # Only need to loop until we actually find an
            # original_file with a successful job.
            if original_file:
                break

            pjs_for_sample = sample.get_processor_jobs()
            for processor_job in list(pjs_for_sample):
                if processor_job.success:
                    original_file = processor_job.original_files.first()
                    if original_file:
                        break

        if not original_file:
            logger.error(
                "Could not find a single sample in the experiment that had a successful Salmon job.",
                experiment=accession_code,
            )
            sys.exit(1)

        job = ProcessorJob()
        job.pipeline_applied = "TXIMPORT"
        job.save()

        pjofa = ProcessorJobOriginalFileAssociation()
        pjofa.processor_job = job
        pjofa.original_file = original_file
        pjofa.save()

        tximport.tximport(job.id)

        sys.exit(0)

import sys
import time

from django.core.management.base import BaseCommand

from data_refinery_common.enums import ProcessorPipeline
from data_refinery_common.job_lookup import determine_ram_amount
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    Dataset,
    Experiment,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
)
from data_refinery_common.utils import queryset_page_iterator

logger = get_and_configure_logger(__name__)


def create_job_for_organism(organism: Organism):
    """Returns a quantpendia job for the provided organism."""
    job = ProcessorJob()
    job.pipeline_applied = ProcessorPipeline.CREATE_QUANTPENDIA.value
    job.save()

    dset = Dataset()
    dset.data = build_dataset(organism)
    dset.scale_by = "NONE"
    dset.aggregate_by = "EXPERIMENT"
    dset.quantile_normalize = False
    dset.quant_sf_only = True
    dset.svd_algorithm = "NONE"
    dset.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = job
    pjda.dataset = dset
    pjda.save()

    # Have to call this after setting the dataset since it's used in
    # the caclulation.
    job.ram_amount = determine_ram_amount(job)
    job.save()

    return job


def build_dataset(organism: Organism):
    data = {}
    experiments = Experiment.objects.filter(
        organisms=organism,
        technology="RNA-SEQ",
    ).distinct()

    for experiment_page in queryset_page_iterator(experiments):
        for experiment in experiment_page:
            # only include the samples from the target organism that have quant.sf files
            experiment_samples = experiment.samples.filter(organism=organism, technology="RNA-SEQ")
            # split the query into two so to avoid timeouts.
            # assume processed rna-seq samples have a quant.sf file
            processed_samples_with_quantsf = experiment_samples.filter(
                is_processed=True
            ).values_list("accession_code", flat=True)
            # and only check for quant file for unprocessed samples
            unprocessed_samples_with_quantsf = (
                experiment_samples.filter(
                    is_processed=False, results__computedfile__filename="quant.sf"
                )
                .values_list("accession_code", flat=True)
                .distinct()
            )

            sample_accession_codes = list(processed_samples_with_quantsf) + list(
                unprocessed_samples_with_quantsf
            )

            if sample_accession_codes:
                data[experiment.accession_code] = sample_accession_codes

        time.sleep(5)

    return data


def create_quantpendia(organisms, organisms_exclude):
    all_organisms = Organism.objects.all()
    if organisms:
        organisms = organisms.upper().replace(" ", "_").split(",")
        all_organisms = all_organisms.filter(name__in=organisms)

    if organisms_exclude:
        organisms = organisms_exclude.upper().replace(" ", "_").split(",")
        all_organisms = all_organisms.exclude(name__in=organisms)

    logger.debug("Generating quantpendia for organisms", organisms=all_organisms)

    created_jobs = []
    for organism in all_organisms:
        # only generate the quantpendia for organisms that have some samples
        # with quant.sf files.
        has_quantsf_files = organism.sample_set.filter(
            technology="RNA-SEQ", results__computedfile__filename="quant.sf"
        ).exists()
        if not has_quantsf_files:
            continue

        job = create_job_for_organism(organism)
        logger.info(
            "Sending compendia job for Organism", job_id=str(job.pk), organism=str(organism)
        )
        send_job(ProcessorPipeline.CREATE_QUANTPENDIA, job)

        created_jobs.append(job)

    return created_jobs


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--organisms", type=str, help=("Comma separated list of organism names.")
        )

        parser.add_argument(
            "--organisms-exclude",
            type=str,
            help=("Comma separated list of organism names that we want to exclude from the list"),
        )

    def handle(self, *args, **options):
        """Create a quantpendia for one or more organisms."""
        create_quantpendia(options["organisms"], options["organisms_exclude"])

        sys.exit(0)

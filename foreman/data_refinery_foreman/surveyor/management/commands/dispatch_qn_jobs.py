from django.core.management.base import BaseCommand
from django.db.models import Count

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    Dataset,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
)

logger = get_and_configure_logger(__name__)

MIN = 100


def dispatch_qn_job_if_eligible(organism: Organism) -> None:
    """Checks if the organism is elgible for a QN job and if so dispatches it.

    An organism is eligible for a QN job if it has more than MIN
    samples on a single platform.
    """
    samples = Sample.processed_objects.filter(
        organism=organism,
        has_raw=True,
        technology="MICROARRAY",
        is_processed=True,
        platform_name__contains="Affymetrix",
    )

    if samples.count() < MIN:
        logger.info(
            "Total proccessed samples don't meet minimum threshhold",
            organism=organism,
            count=samples.count(),
            min=MIN,
        )
        return

    platform_counts = (
        samples.values("platform_accession_code")
        .annotate(dcount=Count("platform_accession_code"))
        .order_by("-dcount")
    )
    biggest_platform = platform_counts[0]["platform_accession_code"]

    sample_codes_results = Sample.processed_objects.filter(
        platform_accession_code=biggest_platform,
        has_raw=True,
        technology="MICROARRAY",
        organism=organism,
        is_processed=True,
    ).values("accession_code")

    if sample_codes_results.count() < MIN:
        logger.info(
            "Number of processed samples for largest platform didn't mean threshold.",
            organism=organism,
            platform_accession_code=biggest_platform,
            count=sample_codes_results.count(),
            min=MIN,
        )
        return

    sample_codes = [res["accession_code"] for res in sample_codes_results]

    dataset = Dataset()
    dataset.data = {organism.name + "_(" + biggest_platform + ")": sample_codes}
    dataset.aggregate_by = "ALL"
    dataset.scale_by = "NONE"
    dataset.quantile_normalize = False
    dataset.save()

    job = ProcessorJob()
    job.pipeline_applied = "QN_REFERENCE"
    job.save()

    pjda = ProcessorJobDatasetAssociation()
    pjda.processor_job = job
    pjda.dataset = dataset
    pjda.save()

    logger.info("Sending QN_REFERENCE for Organism", job_id=str(job.pk), organism=str(organism))
    send_job(ProcessorPipeline.QN_REFERENCE, job)


class Command(BaseCommand):
    def handle(self, *args, **options):
        """Dispatch QN_REFERENCE creation jobs for all Organisms with a
        platform with enough processed samples."""

        if options["organisms"]:
            organism_names = options["organisms"].split(",")
            organisms = Organism.objects.filter(name__in=organism_names)
        else:
            organisms = Organism.objects.all()

        for organism in organisms:
            dispatch_qn_job_if_eligible(organism)

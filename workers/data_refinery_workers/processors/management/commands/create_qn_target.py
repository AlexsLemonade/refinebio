import sys

from django.core.management.base import BaseCommand
from django.db.models import Count

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    Dataset,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
)
from data_refinery_workers.processors import qn_reference

logger = get_and_configure_logger(__name__)


def get_biggest_platform(organism):
    samples = organism.sample_set.filter(has_raw=True, technology="MICROARRAY", is_processed=True)
    if samples.count() == 0:
        logger.error("No processed samples for organism.", organism=organism, count=samples.count())
        return None

    platform_counts = (
        samples.values("platform_accession_code")
        .annotate(dcount=Count("platform_accession_code"))
        .order_by("-dcount")
    )
    return platform_counts[0]["platform_accession_code"]


def create_qn_target(organism, platform, create_results=True):
    sample_codes_results = Sample.processed_objects.filter(
        platform_accession_code=platform,
        has_raw=True,
        technology="MICROARRAY",
        organism=organism,
        is_processed=True,
    ).values("accession_code")
    sample_codes = [res["accession_code"] for res in sample_codes_results]

    dataset = Dataset()
    dataset.data = {organism.name + "_(" + platform + ")": sample_codes}
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

    return qn_reference.create_qn_reference(job.pk, create_results=create_results)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--organism", type=str, help=("Name of organism"))
        parser.add_argument(
            "--all",
            default=False,
            action="store_true",
            help=("Perform for all organisms that meet <min> processed samples"),
        )
        parser.add_argument(
            "--min", type=int, default=100, help=("Minimum number of processed samples")
        )
        parser.add_argument("--platform", type=str, default=None, help=("Name of platform"))
        parser.add_argument("--job_id", type=int, default=None, help=("ID of job to run."))

    def handle(self, *args, **options):
        """ """

        if not options["job_id"]:
            if options["organism"] is None and not options["all"]:
                logger.error("You must specify an organism or --all")
                sys.exit(1)

            if options["organism"] and (options.get("organism", "") != "ALL"):
                organisms = [Organism.get_object_for_name(options["organism"].upper())]
            else:
                organisms = Organism.objects.all()

            for organism in organisms:
                if not organism_can_have_qn_target(organism, options["min"]):
                    logger.error(
                        "Organism does not have any platform with enough samples to generate a qn target",
                        organism=organism,
                        min=options["min"],
                    )
                    continue

                if options["platform"] is None:
                    biggest_platform = get_biggest_platform(organism)
                    if biggest_platform is None:
                        logger.error("No processed samples for organism.", organism=organism)
                        continue
                else:
                    biggest_platform = options["platform"]

                final_context = create_qn_target(organism, platform=biggest_platform)

                if final_context["success"]:
                    print(":D")
                    self.stdout.write("Target file: " + final_context["target_file"])
                    self.stdout.write(
                        "Target S3: " + str(final_context["computed_files"][0].get_s3_url())
                    )
                else:
                    print(":(")
        else:
            qn_reference.create_qn_reference(options["job_id"])


def organism_can_have_qn_target(organism: Organism, sample_threshold=100):
    """Check that the organism has more than `sample_threshold` samples on
    some microarray platform"""
    microarray_platforms = (
        organism.sample_set.filter(has_raw=True, technology="MICROARRAY", is_processed=True)
        .values("platform_accession_code")
        .annotate(count=Count("id"))
        .filter(count__gt=sample_threshold)
    )

    return microarray_platforms.exists()

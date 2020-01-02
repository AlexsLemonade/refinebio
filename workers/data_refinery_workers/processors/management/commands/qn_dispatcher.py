import boto3
import botocore
import nomad
import uuid

from django.core.management.base import BaseCommand
from django.db.models import Count
from nomad.api.exceptions import URLNotFoundNomadException

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_lookup import ProcessorPipeline, SurveyJobTypes
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import SurveyJob, SurveyJobKeyValue
from data_refinery_common.utils import parse_s3_url, get_env_variable

from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SampleComputedFileAssociation,
)

logger = get_and_configure_logger(__name__)

MIN = 100


class Command(BaseCommand):
    def handle(self, *args, **options):
        """ Dispatch QN_REFERENCE creation jobs for all Organisms with a platform with enough processed samples. """

        organisms = Organism.objects.all()

        for organism in organisms:
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
                continue

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
                continue

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

            logger.info(
                "Sending QN_REFERENCE for Organism", job_id=str(job.pk), organism=str(organism)
            )
            send_job(ProcessorPipeline.QN_REFERENCE, job)

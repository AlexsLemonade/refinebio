import time

from django.core.management.base import BaseCommand

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Sample
from data_refinery_foreman.surveyor.sra import SraSurveyor

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--limit",
            default=1000,
            type=int,
            help="Number of samples to refresh",
        )
        parser.add_argument(
            "--source",
            choices=("SRA",),
            required=True,
            type=str,
            help="Source name (ARRAY_EXPRESS, GEO, SRA)",
        )

    def handle(self, *args, **options):
        for sample in Sample.objects.filter(
            developmental_stage__isnull=True,
            last_refreshed__isnull=True,
            source_database=options["source"],
        ).order_by("id")[: options["limit"]]:
            logger.info(f"Refreshing metadata for a sample {sample.accession_code}")
            try:
                _, sample_metadata = SraSurveyor.gather_all_metadata(sample.accession_code)
                SraSurveyor._apply_harmonized_metadata_to_sample(sample_metadata)
            except Exception as e:
                logger.exception(e)
            finally:
                sample.save()

            time.sleep(1)

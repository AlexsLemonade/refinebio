from django.core.management.base import BaseCommand

from data_refinery_common.models import (
    Sample,
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_foreman.surveyor.harmony import harmonize

logger = get_and_configure_logger(__name__)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--source-database",
            type=str,
            help=("The name of a source database, such as Array Express, GEO, or SRA."
                  "All samples from this source database will have their metadata refreshed.")
        )

    def handle(self, *args, **options):
        """Refreshes the metadata for RNA-Seq samples
        """
        if options["source_database"] is None:
            sra_samples = Sample.objects.filter(
                technology="RNA-SEQ"
            )
        else:
            source_database = options["source_database"]
            sra_samples = Sample.objects.filter(
                source_database=source_database
            )

        for sample in sra_samples:
            logger.debug("Refreshing metadata for a sample.",
                         sample=sample.accession_code)
            metadata = SraSurveyor.gather_all_metadata(sample.accession_code)
            SraSurveyor._apply_harmonized_metadata_to_sample(sample, metadata)
            sample.save()

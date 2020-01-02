import sys
import time

from django.core.management.base import BaseCommand

import GEOparse

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Sample
from data_refinery_common.performant_pagination.pagination import PerformantPaginator
from data_refinery_foreman.surveyor import harmony
from data_refinery_foreman.surveyor.array_express import ArrayExpressSurveyor
from data_refinery_foreman.surveyor.geo import GeoSurveyor
from data_refinery_foreman.surveyor.sra import SraSurveyor

logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--source-database",
            type=str,
            help=(
                "The name of a source database, such as ARRAY_EXPRESS, GEO, or SRA."
                "All samples from this source database will have their metadata refreshed."
            ),
        )

    def handle(self, *args, **options):
        """Refreshes the metadata for all samples, or samples from a specific database
        """
        possible_source_databases = ["ARRAY_EXPRESS", "GEO", "SRA"]

        if options.get("source_database", None) is None:
            samples = Sample.objects.all()
        elif options["source_database"] in possible_source_databases:
            source_database = options["source_database"]
            samples = Sample.objects.filter(source_database=source_database)
        else:
            logger.error(
                'Invalid source database "{}"'.format(options["source_database"])
                + "\nPossible source databases: {}".format(", ".join(possible_source_databases))
            )
            sys.exit(1)

        paginator = PerformantPaginator(samples, PAGE_SIZE)
        page = paginator.page()

        while True:
            for sample in samples:
                logger.debug("Refreshing metadata for a sample.", sample=sample.accession_code)
                if sample.source_database == "SRA":
                    metadata = SraSurveyor.gather_all_metadata(sample.accession_code)
                    SraSurveyor._apply_harmonized_metadata_to_sample(sample, metadata)
                elif sample.source_database == "GEO":
                    gse = GEOparse.get_GEO(
                        sample.experiments.first().accession_code,
                        destdir="/tmp/management",
                        how="brief",
                        silent=True,
                    )
                    preprocessed_samples = harmony.preprocess_geo(gse.gsms.items())
                    harmonized_samples = harmony.harmonize(preprocessed_samples)
                    GeoSurveyor._apply_harmonized_metadata_to_sample(
                        sample, harmonized_samples[sample.title]
                    )
                elif sample.source_database == "ARRAY_EXPRESS":
                    SDRF_URL_TEMPLATE = (
                        "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.sdrf.txt"
                    )
                    sdrf_url = SDRF_URL_TEMPLATE.format(
                        code=sample.experiments.first().accession_code
                    )
                    sdrf_samples = harmony.parse_sdrf(sdrf_url)
                    harmonized_samples = harmony.harmonize(sdrf_samples)
                    ArrayExpressSurveyor._apply_harmonized_metadata_to_sample(
                        sample, harmonized_samples[sample.title]
                    )

                sample.save()

            if not page.has_next():
                break
            else:
                page = paginator.page(page.next_page_number())

            # 2000 samples queued up every five minutes should be fast
            # enough and also not thrash the DB.
            time.sleep(60 * 5)

"""
This command checks if samples from GEO have incorrect platform
information and re-surveys them.

We used to be worse at surveying from GEO so some samples have bad
metadata and some of them were even successufully processed. Therefore
if this command find any with incorrect metadata it will purge the
experiment and all related data from the database so we can re-survey
the experiment and reprocess it correctly.
"""

import GEOparse
import os

from django.core.management.base import BaseCommand
from django.utils import timezone

from data_refinery_common.models import Experiment, Sample, CdfCorrectedAccession
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_internal_microarray_accession
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment


logger = get_and_configure_logger(__name__)

PAGE_SIZE = 2000


GEO_TEMP_DIR = "/tmp/"


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            help="Prints which experiments would be resurveyed without doing so.",
            action="store_true",
        )

    def handle(self, *args, **options):
        """Re-surveys GEO experiments containing samples with incorrect platform information.
        """
        # Check against CDF corrected accessions table to prevent recorrection of the same samples.
        corrected_experiments = CdfCorrectedAccession.objects.all().values("accession_code")
        corrected_accessions = [
            experiment["accession_code"] for experiment in corrected_experiments
        ]

        gse_experiments = Experiment.objects.filter(source_database="GEO").exclude(
            accession_code__in=corrected_experiments
        )

        paginator = Paginator(gse_experiments, PAGE_SIZE)
        page = paginator.page()

        while True:
            for experiment in page.object_list:
                try:
                    gse = GEOparse.get_GEO(
                        experiment.accession_code, destdir=GEO_TEMP_DIR, how="brief", silent=True
                    )

                    sample_accessions = list(gse.gsms.keys())
                    samples = Sample.objects.filter(accession_code__in=sample_accessions)

                    wrong_platform = False
                    for sample in samples:
                        gpl = gse.gsms[sample.accession_code].metadata["platform_id"][0]
                        internal_accession = get_internal_microarray_accession(gpl)
                        if internal_accession != sample.platform_accession_code:
                            wrong_platform = True
                            break

                    if wrong_platform:
                        if options["dry_run"]:
                            logger.info(
                                "Would have re-surveyed experiment with accession code %s",
                                experiment.accession_code,
                            )
                        else:
                            logger.info(
                                "Re-surveying experiment with accession code %s",
                                experiment.accession_code,
                            )

                            purge_experiment(experiment.accession_code)

                            queue_surveyor_for_accession(experiment.accession_code)

                    current_time = timezone.now()
                    CdfCorrectedAccession(
                        accession_code=experiment.accession_code, created_at=current_time
                    ).save()
                except Exception:
                    logger.exception("Caught an exception with %s!", experiment.accession_code)
                finally:
                    # GEOparse downloads files here and never cleans them up! Grrrr!
                    download_path = GEO_TEMP_DIR + experiment.accession_code + "_family.soft.gz"
                    # It's not a directory, but ignore_errors is useful.
                    try:
                        os.remove(download_path)
                    except:
                        # Don't anything interrupt this, like say,
                        # GEOParse downloading a directory instead of
                        # a file...
                        logger.exception("Failed to delete an archive.")

            if not page.has_next():
                break

            page = paginator.page(page.next_page_number())

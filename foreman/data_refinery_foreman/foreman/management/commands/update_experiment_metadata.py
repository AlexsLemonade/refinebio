import sys
import time

from django.core.management.base import BaseCommand

import GEOparse

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Experiment
from data_refinery_common.performant_pagination.pagination import PerformantPaginator
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.array_express import EXPERIMENTS_URL, ArrayExpressSurveyor
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
                "All experiments from this source database will have their metadata refreshed."
            ),
        )

    def handle(self, *args, **options):
        """Refreshes the metadata for all experiments, or experiments from a specific database
        """
        possible_source_databases = ["ARRAY_EXPRESS", "GEO", "SRA"]

        if options.get("source_database", None) is None:
            experiments = Experiment.objects.all()
        elif options["source_database"] in possible_source_databases:
            source_database = options["source_database"]
            experiments = Experiment.objects.filter(source_database=source_database)
        else:
            logger.error(
                'Invalid source database "{}"'.format(options["source_database"])
                + "\nPossible source databases: {}".format(", ".join(possible_source_databases))
            )
            sys.exit(1)

        paginator = PerformantPaginator(experiments, PAGE_SIZE)
        page = paginator.page()

        while True:
            for experiment in page.object_list:
                logger.debug(
                    "Refreshing metadata for an experiment.", experiment=experiment.accession_code
                )
                try:
                    if experiment.source_database == "SRA":
                        metadata = SraSurveyor.gather_all_metadata(
                            experiment.samples.first().accession_code
                        )
                        SraSurveyor._apply_metadata_to_experiment(experiment, metadata)

                    elif experiment.source_database == "GEO":
                        gse = GEOparse.get_GEO(
                            experiment.accession_code,
                            destdir="/tmp/management",
                            how="brief",
                            silent=True,
                        )

                        GeoSurveyor._apply_metadata_to_experiment(experiment, gse)

                    elif experiment.source_database == "ARRAY_EXPRESS":
                        request_url = EXPERIMENTS_URL + experiment.accession_code
                        experiment_request = utils.requests_retry_session().get(
                            request_url, timeout=60
                        )
                        try:
                            parsed_json = experiment_request.json()["experiments"]["experiment"][0]
                        except KeyError:
                            logger.error(
                                "Remote experiment has no Experiment data!",
                                experiment_accession_code=experiment.accession_code,
                                survey_job=self.survey_job.id,
                            )
                            continue
                        ArrayExpressSurveyor._apply_metadata_to_experiment(experiment, parsed_json)

                    experiment.save()

                # If there are any errors, just continue. It's likely that it's
                # just a problem with this experiment.
                except Exception:
                    logger.exception(
                        "exception caught while updating metadata for {}".format(
                            experiment.accession_code
                        )
                    )

            if not page.has_next():
                break
            else:
                page = paginator.page(page.next_page_number())

            # 2000 samples queued up every five minutes should be fast
            # enough and also not thrash the DB.
            time.sleep(60 * 5)

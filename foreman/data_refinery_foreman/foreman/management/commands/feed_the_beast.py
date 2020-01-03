"""This command will slowly retry Salmon jobs that timed out.
This is now necessary because samples with unmated reads will no longer cause
us to time out. It will only queue 300 an hour so as to not overload ENA.
"""

import time
from typing import List

from django.core.management.base import BaseCommand
from django.utils import timezone
from nomad import Nomad

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Experiment, ProcessorJob, SurveyedAccession
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator
from data_refinery_common.utils import get_env_variable
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def handle(self, *args, **options):
        nomad_host = get_env_variable("NOMAD_HOST")
        nomad_port = get_env_variable("NOMAD_PORT", "4646")
        nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)

        with open("config/all_rna_seq_accessions.txt") as accession_list_file:
            all_rna_accessions = [line.strip() for line in accession_list_file]

        with open("config/all_microarray_accessions.txt") as accession_list_file:
            all_microarray_accessions = [line.strip() for line in accession_list_file]

        all_accessions = all_microarray_accessions + all_rna_accessions

        BATCH_SIZE = 1000
        batch_index = 0
        batch_accessions = all_accessions[0:BATCH_SIZE]

        fed_accessions = []

        while batch_accessions:
            logger.info(
                "Looping through another batch of 1000 experiments, starting with accession code: %s",
                batch_accessions[0],
            )

            # Check against surveyed accessions table to prevent resurveying
            surveyed_experiments = SurveyedAccession.objects.filter(
                accession_code__in=batch_accessions
            ).values("accession_code")

            surveyed_accessions = [
                experiment["accession_code"] for experiment in surveyed_experiments
            ]

            missing_accessions = set(batch_accessions) - set(surveyed_accessions)
            while len(missing_accessions) > 0:
                try:
                    all_surveyor_jobs = nomad_client.jobs.get_jobs(prefix="SURVEYOR")

                    num_surveyor_jobs = 0
                    for job in all_surveyor_jobs:
                        if job["ParameterizedJob"] and job["JobSummary"].get("Children", None):
                            num_surveyor_jobs = (
                                num_surveyor_jobs + job["JobSummary"]["Children"]["Pending"]
                            )
                            num_surveyor_jobs = (
                                num_surveyor_jobs + job["JobSummary"]["Children"]["Running"]
                            )
                except:
                    logger.exception("Exception caught counting surveyor jobs!")
                    # Probably having trouble communicating with Nomad, let's try again next loop.
                    continue

                if num_surveyor_jobs < 15:
                    accession_code = missing_accessions.pop()
                    try:
                        queue_surveyor_for_accession(accession_code)
                        fed_accessions.append(accession_code)
                        time.sleep(30)
                    except:
                        # We don't want to stop, gotta keep feeding the beast!!!!
                        logger.exception(
                            "Exception caught while looping through all accessions!",
                            accession_code=accession_code,
                        )
                else:
                    # Do it here so we don't sleep when there's an exception
                    time.sleep(30)

            # Bulk insert fed_accessions to SurveyedAccession
            new_surveyed_accessions = []
            current_time = timezone.now()

            for accession in fed_accessions:
                new_surveyed_accessions.append(
                    SurveyedAccession(accession_code=accession, created_at=current_time)
                )

            SurveyedAccession.objects.bulk_create(new_surveyed_accessions)
            fed_accessions = []

            batch_index += 1
            if batch_index * BATCH_SIZE >= len(all_accessions):
                break

            batch_start = batch_index * BATCH_SIZE
            batch_end = batch_start + BATCH_SIZE
            batch_accessions = all_accessions[batch_start:batch_end]

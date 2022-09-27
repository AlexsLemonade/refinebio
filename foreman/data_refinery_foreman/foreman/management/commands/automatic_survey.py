"""
Automatic survey command.
"""

from django.core.management.base import BaseCommand
from django.template.defaultfilters import pluralize
from django.utils import timezone

from dateutil.relativedelta import relativedelta

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import GatheredAccession, SurveyedAccession, SurveyJob
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("-d", "--days", help="The time window in days.", type=int, default=1)
        parser.add_argument(
            "-q",
            "--quota",
            help="The number of accessions to survey per the time window.",
            type=int,
            default=10,
        )

    def handle(self, *args, **options):
        create_survey_jobs(options["days"], options["quota"])


def create_survey_jobs(days, quota):
    """Automatically creates survey jobs based on the passed time window requirements."""

    remaining_quota = get_remaining_quota(days, quota)
    if not remaining_quota:
        return

    gathered_accessions = [
        ga["accession_code"]
        for ga in GatheredAccession.objects.order_by("created_at").values("accession_code")
    ]

    batch_start = 0
    batch_size = 1000
    queued_accessions = list()
    while True:
        batch_accessions = gathered_accessions[batch_start : batch_start + batch_size]
        if not batch_accessions:
            break

        surveyed_accessions = (
            sa["accession_code"]
            for sa in SurveyedAccession.objects.filter(accession_code__in=batch_accessions).values(
                "accession_code"
            )
        )

        for accession_code in set(batch_accessions) - set(surveyed_accessions):
            try:
                queue_surveyor_for_accession(accession_code)
            except Exception as e:
                logger.error(f"Couldn't queue accession {accession_code} due to: {e}")
            else:
                queued_accessions.append(accession_code)
                if len(queued_accessions) == remaining_quota:
                    break

        if len(queued_accessions) == remaining_quota:
            break

        batch_start += batch_size

    queued_accession_count = len(queued_accessions)
    logger.info(
        f"Queued {queued_accession_count} accession{pluralize(queued_accession_count)} "
        f"based on {quota} accession{pluralize(quota)} per {days} "
        f"day{pluralize(days)} quota."
    )
    SurveyedAccession.objects.bulk_create(
        (SurveyedAccession(accession_code=accession_code) for accession_code in queued_accessions)
    )


def get_remaining_quota(days, quota):
    """Returns the remaining quota or 0 if nothing left."""
    time_window_survey_jobs = SurveyJob.objects.filter(
        created_at__gt=timezone.now() - relativedelta(days=days)
    )

    return max(0, quota - time_window_survey_jobs.count())

# import csv
# import json
import os

from django.conf import settings
from django.core.management.base import BaseCommand

import pandas as pd
import requests
import scipy.stats

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Organism
from data_refinery_workers.processors import qn_reference
from data_refinery_workers.processors.management.commands import create_qn_target

logger = get_and_configure_logger(__name__)

KS_TEST_THRESHOLD = 0.05


def signal_failure(bad_organisms):
    if not settings.RUNNING_IN_CLOUD:
        logger.error("Cannot message Slack when we are not running in the cloud")
        return
    elif not settings.ENGAGEMENTBOT_WEBHOOK:
        logger.error("Cannot message slack because we don't have a valid webhook")
        return

    requests.post(
        settings.ENGAGEMENTBOT_WEBHOOK,
        json={
            "attachments": [
                {
                    "fallback": "QN Reference Test Failed",
                    "title": "QN Reference Test Failed",
                    "color": "#db3b28",
                    "text": "Some newly-generated QN references were not similar-enough to our references currently in use.",
                    "fields": [
                        {
                            "title": "Failing Organisms",
                            "value": ", ".join(o.name for o in bad_organisms),
                        }
                    ],
                    "footer": "Refine.bio",
                    "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                }
            ],
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )


def check_qn_reference_agreement():
    bad_organisms = []
    for organism in Organism.objects.filter(qn_target__isnull=False):
        final_context = create_qn_target.create_qn_target(
            organism,
            platform=create_qn_target.get_biggest_platform(organism),
            create_results=False,
        )
        new_target = pd.read_csv(
            final_context["target_file"],
            index=False,
            index_col=0,
            header=False,
            sep="\t",
            encoding="utf-8",
        )

        existing_target_computedfile = organism.qn_target.computedfile_set.latest()
        existing_target_filename = existing_target_computedfile.sync_from_s3()
        existing_target = pd.read_csv(
            existing_target_filename,
            index=False,
            index_col=0,
            header=False,
            sep="\t",
            encoding="utf-8",
        )
        # We don't want this file hanging around on the foreman for too long
        existing_target_computedfile.delete_local_file()

        _, pvalue = scipy.stats.ks_2samp(new_target, existing_target)

        # TODO: decide on a good threshold
        if pvalue < 0.05:
            bad_organisms.append(organism)

    if len(bad_organisms) != 0:
        signal_failure(bad_organisms)

        pass


class Command(BaseCommand):
    def handle(self, *args, **options):
        try:
            check_qn_reference_agreement()
        except Exception as e:
            requests.post(
                settings.ENGAGEMENTBOT_WEBHOOK,
                json={
                    "attachments": [
                        {
                            "fallback": "Exception raised during QN Reference Test",
                            "title": "Exception raised during QN Reference Test",
                            "color": "#db3b28",
                            "text": f"```\n{str(e)}\n```",
                            "footer": "Refine.bio",
                            "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                        }
                    ],
                },
                headers={"Content-Type": "application/json"},
                timeout=10,
            )

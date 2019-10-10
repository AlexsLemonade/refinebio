import random
import sys
import time
import re
from typing import Dict, List

from django.core.management.base import BaseCommand
from django.db.models import OuterRef, Subquery, Count
from dateutil.parser import parse as parse_date
from data_refinery_common.models import (
    ProcessorJob,
    Sample,
    ComputedFile
)
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.job_management import create_downloader_job
from data_refinery_common.performant_pagination.pagination import PerformantPaginator as Paginator


logger = get_and_configure_logger(__name__)

def remove_qn_targets():
    pass

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            '--qn-targets',
            action='store_true',
            help='Finds any invalid QN Targets that got generated and removes them from the db.',
        )

    def handle(self, *args, **options):
        """ Re-queues all unprocessed RNA-Seq samples for an organism. """
        if options["qn_targets"]:
            remove_qn_targets()


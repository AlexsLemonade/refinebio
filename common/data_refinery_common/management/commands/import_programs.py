#!/usr/bin/env python

"""Read an input tab-separated file and add each line into Program model.
This management commands requires an input filename, and the envioronment
variable $CIRCLE_TAG must be defined with its value starting with "v".
Each line in <filename> must have exactly three columns delimited by tab.
All lines that are blank or start with "#" will be skipped.

Usage:
  python manage.py import_programs <filename>

To use your own $CIRCLE_TAG, launch the command like this:
  CIRCLE_TAG=<my_tag> python manage.py import_programs <filename>
"""

import argparse
import os
from django.core.management.base import BaseCommand
from django.db import transaction
from data_refinery_common.models import Program
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('filename', type=argparse.FileType('r'))
        help = ("Input tab-separated file that lists programs.")

    def handle(self, **options):
        # Verify that $CIRCLE_TAG a valid envioronment variable
        try:
            git_tag = os.environ['CIRCLE_TAG']
        except:
            logger.error("CIRCLE_TAG not defined")
            return 1
        # Verify that $CIRCLE_TAG starts with 'v'.
        if not git_tag.startswith('v'):
            logger.error("CIRCLE_TAG must start with 'v'")
            return 1

        try:
            import_data(options['filename'], git_tag)
            logger.info(self.style.NOTICE("Program data imported successfully"))
        except Exception as e:
            logger.error("import_programs raised an exception: %s", e)
            return 1

        return 0


def import_data(file_handle, git_tag):
    """Read each line from file_handle into "programs" table.
    The input git_tag will be "system_version" in the table.
    """

    with transaction.atomic():
        for line_index, line in enumerate(file_handle):
            line = line.strip()
              # Skip lines that are blank or start with "#"
            if len(line) == 0 or line.startswith('#'):
                continue

            tokens = line.split('\t')
            if len(tokens) != 3:
                raise Exception("Input file line #%: number of columns is not 3" %
                                (line_index + 1))

            Program.objects.create(name=tokens[0],
                                   version=tokens[1],
                                   command=tokens[2],
                                   system_version=git_tag)

import os
import psutil

from django.test import TestCase, tag
from typing import List
from unittest.mock import patch, call
from urllib.error import URLError

from data_refinery_workers.downloaders import utils

class UtilsTestCase(TestCase):
    @tag('downloaders')
    def test_no_jobs_to_create(self):
        """Make sure this function doesn't raise an exception with no files."""
        utils.create_processor_job_for_original_files([])

        self.assertTrue(True)

    @tag('downloaders')
    def test_max_jobs(self):
        """ Test max jobs before job suicide """

        total_vm = psutil.virtual_memory().total
        gb = int(total_vm / 1000000000)

        max_jobs = utils.get_max_jobs_for_current_node()
        # We're not going to run our tests on a prod box, so this should always be True.
        self.assertNotEqual(max_jobs, 8)
        self.assertNotEqual(max_jobs, None)

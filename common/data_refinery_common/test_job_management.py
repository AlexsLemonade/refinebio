
from django.test import TestCase, tag
from typing import List
from unittest.mock import patch, call
from urllib.error import URLError

from data_refinery_common.job_management import create_processor_job_for_original_files


class UtilsTestCase(TestCase):
    def test_no_jobs_to_create(self):
        """Make sure this function doesn't raise an exception with no files."""
        create_processor_job_for_original_files([])

        self.assertTrue(True)

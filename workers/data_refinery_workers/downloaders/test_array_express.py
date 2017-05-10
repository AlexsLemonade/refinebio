from django.test import TestCase
from data_refinery_models.models import Batch
from data_refinery_workers.downloaders import array_express


class DownloadArrayExpressTestCase(TestCase):
    def test_good_batch_grouping(self):
        """Returns true if all batches have the same download_url."""
        batches = [Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com")]
        job_id = 1
        self.assertIsNone(array_express._verify_batch_grouping(batches, job_id))

    def test_bad_batch_grouping(self):
        """Returns true if all batches have the same download_url."""
        batches = [Batch(download_url="https://example.com"),
                   Batch(download_url="https://example.com"),
                   Batch(download_url="https://wompwomp.com")]
        job_id = 1
        with self.assertRaises(ValueError):
            array_express._verify_batch_grouping(batches, job_id)

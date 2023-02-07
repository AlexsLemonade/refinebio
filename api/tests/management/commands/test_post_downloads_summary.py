from datetime import timedelta
from unittest.mock import patch

from django.test import TestCase
from django.utils import timezone

from data_refinery_api.management.commands.post_downloads_summary import post_downloads_summary
from data_refinery_common.models import Dataset, DatasetAnnotation


class DownloadsPostTestCase(TestCase):
    @patch("data_refinery_common.models.organism.requests.post")
    def test_post_download_summary_one_download(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville"}).save()

        post_downloads_summary(7, "ccdl-general")

        text = requests_post.mock_calls[0][2]["json"]["text"]
        self.assertIn("1 users downloaded 1 datasets from 1 locations.", text)

    @patch("data_refinery_common.models.organism.requests.post")
    def test_post_download_summary_two_downloads(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville"}).save()
        dataset = Dataset(email_address="test2@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Not Testyville"}).save()

        post_downloads_summary(7, "ccdl-general")

        text = requests_post.mock_calls[0][2]["json"]["text"]
        self.assertIn("2 users downloaded 2 datasets from 2 locations.", text)

    @patch("data_refinery_common.models.organism.requests.post")
    def test_post_download_summary_filtered(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville"}).save()
        last_week = timezone.now() - timedelta(days=8)
        dataset = Dataset(created_at=last_week, email_address="test2@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(
            created_at=last_week, dataset=dataset, data={"location": "Not Testyville"}
        ).save()

        post_downloads_summary(7, "ccdl-general")

        text = requests_post.mock_calls[0][2]["json"]["text"]
        self.assertIn("1 users downloaded 1 datasets from 1 locations.", text)

    @patch("data_refinery_common.models.organism.requests.post")
    def test_empty_post_download_summary(self, requests_post):
        post_downloads_summary(7, "ccdl-general")

        text = requests_post.mock_calls[0][2]["json"]["text"]
        self.assertIn("no downloads", text)

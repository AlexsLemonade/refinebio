from datetime import timedelta
from unittest.mock import patch

from django.test import TestCase
from django.utils import timezone

from data_refinery_api.management.commands.post_downloads_summary import post_downloads_summary
from data_refinery_common.models import Dataset, DatasetAnnotation


class DownloadsPostTestCase(TestCase):
    @patch("data_refinery_api.management.commands.post_downloads_summary.requests.post")
    def test_post_download_summary_one_download(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville, TSA"}).save()

        post_downloads_summary(7, "ccdl-general-test")
        json = requests_post.mock_calls[0][2]["json"]

        self.assertEqual(
            "In the last 7 days, 1 user downloaded 1 dataset from 1 location.", json["text"]
        )
        self.assertEqual(
            "\n".join(("*New users*", "test@gmail.com | 1 download from Testyville, TSA")),
            json["blocks"][1]["text"]["text"],
        )
        self.assertEqual(
            "\n".join(("*Top 1 country*", "TSA: 1 download")), json["blocks"][2]["text"]["text"]
        )

    @patch("data_refinery_api.management.commands.post_downloads_summary.requests.post")
    def test_post_download_summary_two_downloads(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville, TSA"}).save()

        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville, TSA"}).save()

        dataset = Dataset(email_address="test2@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Not Testyville, Not TSA"}).save()

        post_downloads_summary(7, "ccdl-general-test")
        json = requests_post.mock_calls[0][2]["json"]

        self.assertEqual(
            "In the last 7 days, 2 users downloaded 3 datasets from 2 locations.", json["text"]
        )
        self.assertEqual(
            "\n".join(
                (
                    "*New users*",
                    "test@gmail.com | 2 downloads from Testyville, TSA",
                    "test2@gmail.com | 1 download from Not Testyville, Not TSA",
                )
            ),
            json["blocks"][1]["text"]["text"],
        )
        self.assertEqual(
            "\n".join(("*Top 2 countries*", "TSA: 2 downloads", "Not TSA: 1 download")),
            json["blocks"][2]["text"]["text"],
        )

    @patch("data_refinery_api.management.commands.post_downloads_summary.requests.post")
    def test_post_download_summary_ordering(self, requests_post):
        for _ in range(2):
            dataset = Dataset(email_address="ATEST2@gmail.com", is_processed=True)
            dataset.save()
            DatasetAnnotation(dataset=dataset, data={"location": "Not Testyville, Not TSA"}).save()

        for _ in range(2):
            dataset = Dataset(email_address="TEST@gmail.com", is_processed=True)
            dataset.save()
            DatasetAnnotation(dataset=dataset, data={"location": "Testyville, TSA"}).save()

        dataset = Dataset(email_address="TEST3@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testbec City, Tanada"}).save()

        for _ in range(3):
            dataset = Dataset(email_address="WEST@gmail.com", is_processed=True)
            dataset.save()
            DatasetAnnotation(dataset=dataset, data={"location": "Westyville, Texico"}).save()

        post_downloads_summary(7, "ccdl-general-test")
        json = requests_post.mock_calls[0][2]["json"]

        self.assertEqual(
            "In the last 7 days, 4 users downloaded 8 datasets from 4 locations.", json["text"]
        )
        self.assertEqual(
            "\n".join(
                (
                    "*New users*",
                    "west@gmail.com | 3 downloads from Westyville, Texico",
                    "atest2@gmail.com | 2 downloads from Not Testyville, Not TSA",
                    "test@gmail.com | 2 downloads from Testyville, TSA",
                    "test3@gmail.com | 1 download from Testbec City, Tanada",
                )
            ),
            json["blocks"][1]["text"]["text"],
        )
        self.assertEqual(
            "\n".join(
                (
                    "*Top 4 countries*",
                    "Texico: 3 downloads",
                    "Not TSA: 2 downloads",
                    "TSA: 2 downloads",
                    "Tanada: 1 download",
                )
            ),
            json["blocks"][2]["text"]["text"],
        )

    @patch("data_refinery_api.management.commands.post_downloads_summary.requests.post")
    def test_post_download_summary_filtered(self, requests_post):
        dataset = Dataset(email_address="test@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(dataset=dataset, data={"location": "Testyville, TSA"}).save()

        last_week = timezone.now() - timedelta(days=8)
        dataset = Dataset(created_at=last_week, email_address="test2@gmail.com", is_processed=True)
        dataset.save()
        DatasetAnnotation(
            created_at=last_week, dataset=dataset, data={"location": "Testyville, TSA"}
        ).save()

        post_downloads_summary(7, "ccdl-general-test")
        json = requests_post.mock_calls[0][2]["json"]
        self.assertEqual(
            "In the last 7 days, 1 user downloaded 1 dataset from 1 location.", json["text"]
        )

    @patch("data_refinery_api.management.commands.post_downloads_summary.requests.post")
    def test_empty_post_download_summary(self, requests_post):
        post_downloads_summary(7, "ccdl-general-test")

        json = requests_post.mock_calls[0][2]["json"]
        self.assertEqual("There were no downloads in the last 7 days.", json["text"])

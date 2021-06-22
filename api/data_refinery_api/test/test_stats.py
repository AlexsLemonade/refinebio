from django.urls import reverse
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION


class StatsTestCases(APITestCase):
    """Basic sanity test. Should be improved."""

    def test_stats_empty(self):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)

    def test_stats(self):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)

from unittest.mock import Mock, patch

from django.http import HttpResponseForbidden, HttpResponseServerError
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION
from data_refinery_api.views import DatasetView, ExperimentListView


class APITestCases(APITestCase):
    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_ok(self, mock_client):
        # We don't even import raven if it's a good response.
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        mock_client.is_enabled.assert_not_called()

    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_404(self, mock_client):
        # We don't send anything to raven if it's not enabled
        mock_client.is_enabled.side_effect = lambda: False
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(
            reverse(
                "experiments_detail",
                kwargs={"accession_code": "INEXISTENT", "version": API_VERSION},
            )[:-1]
            + "aasdas/"
        )
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

    @patch.object(ExperimentListView, "get")
    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_403(self, mock_client, mock_get_method):
        mock_get_method.side_effect = Mock(return_value=HttpResponseForbidden())
        # A 403 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = Mock(return_value=True)
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 403)
        mock_client.captureMessage.assert_called()

    @patch.object(ExperimentListView, "get")
    @patch("raven.contrib.django.models.client")
    def test_sentry_middleware_500(self, mock_client, mock_get_method):
        def raise_error(_):
            raise KeyError()

        mock_get_method.side_effect = Mock(return_value=HttpResponseServerError())
        # A 500 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = Mock(return_value=True)
        response = self.client.get(reverse("experiments", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 500)
        mock_client.captureMessage.assert_called()

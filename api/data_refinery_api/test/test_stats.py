from unittest.mock import patch

from django.urls import reverse
from rest_framework.test import APITestCase

from data_refinery_api.test.test_api_general import API_VERSION


MOCK_NOMAD_RESPONSE = [
    {
        "CreateIndex": 5145,
        "ID": "TXIMPORT",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 0},
            "CreateIndex": 5145,
            "JobID": "TXIMPORT",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "TXIMPORT",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_1_2323",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 1},
            "CreateIndex": 5145,
            "JobID": "SALMON_1_2323",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_1_2323",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_3_2534/dispatch-213123-123123-123123",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 0},
            "CreateIndex": 5145,
            "JobID": "SALMON_3_2534/dispatch-213123-123123-123123",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_3_2534/dispatch-213123-123123-123123",
        "ParameterizedJob": True,
        "ParentID": "SALMON_1_2323",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
    {
        "CreateIndex": 5145,
        "ID": "SALMON_2_2323",
        "JobModifyIndex": 5145,
        "JobSummary": {
            "Children": {"Dead": 0, "Pending": 0, "Running": 1},
            "CreateIndex": 5145,
            "JobID": "SALMON_1_2323",
            "ModifyIndex": 5145,
            "Namespace": "default",
            "Summary": {},
        },
        "ModifyIndex": 5145,
        "Name": "SALMON_1_2323",
        "ParameterizedJob": True,
        "ParentID": "",
        "Periodic": False,
        "Priority": 50,
        "Status": "running",
        "StatusDescription": "",
        "Stop": False,
        "SubmitTime": 1552322030836469355,
        "Type": "batch",
    },
]


class StatsTestCases(APITestCase):
    @patch("data_refinery_common.utils.get_nomad_jobs", return_value=[])
    def test_nomad_stats_empty(self, mock_get_nomad_jobs):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["nomad_running_jobs"], 0)
        self.assertEqual(response.json()["nomad_pending_jobs"], 0)

    @patch("data_refinery_common.utils.get_nomad_jobs", return_value=MOCK_NOMAD_RESPONSE)
    def test_nomad_stats(self, mock_get_nomad_jobs):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["nomad_running_jobs"], 2)
        self.assertEqual(response.json()["nomad_pending_jobs"], 0)
        self.assertEqual(response.json()["nomad_running_jobs_by_type"]["SALMON"], 2)
        self.assertEqual(response.json()["nomad_running_jobs_by_volume"]["1"], 1)
        self.assertEqual(response.json()["nomad_running_jobs_by_volume"]["2"], 1)
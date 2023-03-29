from unittest.mock import patch

from django.test import override_settings
from django.urls import reverse
from rest_framework.test import APITestCase

from tests.views.test_api_general import API_VERSION

from data_refinery_api.views.stats import get_batch_jobs_breakdown

QUEUE_NAMES = [
    "data-refinery-batch-compendia-queue-tests-dev",
    "data-refinery-batch-smasher-queue-tests-dev",
    "data-refinery-batch-workers-queue-tests-dev-0",
]


def dummy_get_jobs_in_queue(queue):
    if queue not in QUEUE_NAMES:
        raise ValueError(f"Tried to get jobs for unrecognzied job queue {queue}")
    return {
        # The queues are defined at the bottom of the file because they're pretty long
        "data-refinery-batch-compendia-queue-tests-dev": COMPENDIA_QUEUE,
        "data-refinery-batch-smasher-queue-tests-dev": SMASHER_QUEUE,
        "data-refinery-batch-workers-queue-tests-dev-0": WORKER_QUEUE,
    }[queue]


class StatsTestCases(APITestCase):
    def test_stats_empty(self):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)

    def test_stats(self):
        response = self.client.get(reverse("stats", kwargs={"version": API_VERSION}))
        self.assertEqual(response.status_code, 200)

    @patch("data_refinery_api.views.stats.get_jobs_in_queue")
    @override_settings(AWS_BATCH_QUEUE_ALL_NAMES=QUEUE_NAMES)
    def test_stats_get_batch_breakdown(self, mock_get_jobs_in_queue):
        """Make sure that the batch breakdown has the right stats"""
        mock_get_jobs_in_queue.side_effect = dummy_get_jobs_in_queue

        breakdown = get_batch_jobs_breakdown(force=True)

        self.assertEqual(
            set(breakdown.keys()),
            {
                "pending_jobs",
                "running_jobs",
                "pending_jobs_by_type",
                "running_jobs_by_type",
                "pending_jobs_by_queue",
                "running_jobs_by_queue",
            },
        )

        self.assertEqual(
            breakdown["pending_jobs"],
            PENDING_COMPENDIA_JOBS
            + PENDING_SMASHER_JOBS
            + PENDING_DOWNLOADER_JOBS
            + PENDING_SALMON_JOBS
            + PENDING_AFFY_JOBS,
        )

        self.assertEqual(
            breakdown["running_jobs"],
            RUNNING_COMPENDIA_JOBS
            + RUNNING_SMASHER_JOBS
            + RUNNING_DOWNLOADER_JOBS
            + RUNNING_SALMON_JOBS
            + RUNNING_AFFY_JOBS,
        )

        self.assertEqual(
            breakdown["pending_jobs_by_type"],
            {
                "CREATE_COMPENDIA": PENDING_COMPENDIA_JOBS,
                "SMASHER": PENDING_SMASHER_JOBS,
                "DOWNLOADER": PENDING_DOWNLOADER_JOBS,
                "SALMON": PENDING_SALMON_JOBS,
                "AFFY_TO_PCL": PENDING_AFFY_JOBS,
            },
        )

        self.assertEqual(
            breakdown["running_jobs_by_type"],
            {
                "CREATE_COMPENDIA": RUNNING_COMPENDIA_JOBS,
                "SMASHER": RUNNING_SMASHER_JOBS,
                "DOWNLOADER": RUNNING_DOWNLOADER_JOBS,
                "SALMON": RUNNING_SALMON_JOBS,
                "AFFY_TO_PCL": RUNNING_AFFY_JOBS,
            },
        )

        self.assertEqual(
            breakdown["pending_jobs_by_queue"],
            {
                "data-refinery-batch-compendia-queue-tests-dev": PENDING_COMPENDIA_JOBS,
                "data-refinery-batch-smasher-queue-tests-dev": PENDING_SMASHER_JOBS,
                "data-refinery-batch-workers-queue-tests-dev-0": PENDING_DOWNLOADER_JOBS
                + PENDING_SALMON_JOBS
                + PENDING_AFFY_JOBS,
            },
        )

        self.assertEqual(
            breakdown["running_jobs_by_queue"],
            {
                "data-refinery-batch-compendia-queue-tests-dev": RUNNING_COMPENDIA_JOBS,
                "data-refinery-batch-smasher-queue-tests-dev": RUNNING_SMASHER_JOBS,
                "data-refinery-batch-workers-queue-tests-dev-0": RUNNING_DOWNLOADER_JOBS
                + RUNNING_SALMON_JOBS
                + RUNNING_AFFY_JOBS,
            },
        )


PENDING_COMPENDIA_JOBS = 10
RUNNING_COMPENDIA_JOBS = 1
COMPENDIA_QUEUE = [
    *[
        {"jobName": f"tests_dev_CREATE_COMPENDIA_{i}", "status": "PENDING"}
        for i in range(PENDING_COMPENDIA_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_CREATE_COMPENDIA_{i}", "status": "RUNNING"}
        for i in range(RUNNING_COMPENDIA_JOBS)
    ],
]

PENDING_SMASHER_JOBS = 27
RUNNING_SMASHER_JOBS = 5
# Create some finished jobs that should get ignored
FINISHED_SMASHER_JOBS = 10
SMASHER_QUEUE = [
    *[
        {"jobName": f"tests_dev_SMASHER_{i}", "status": "RUNNABLE"}
        for i in range(PENDING_SMASHER_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_SMASHER_{i}", "status": "RUNNING"}
        for i in range(RUNNING_SMASHER_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_SMASHER_{i}", "status": "SUCCEEDED"}
        for i in range(FINISHED_SMASHER_JOBS)
    ],
]


PENDING_DOWNLOADER_JOBS = 14
RUNNING_DOWNLOADER_JOBS = 10
PENDING_SALMON_JOBS = 2
RUNNING_SALMON_JOBS = 8
PENDING_AFFY_JOBS = 9
RUNNING_AFFY_JOBS = 1
WORKER_QUEUE = [
    *[
        {"jobName": f"tests_dev_DOWNLOADER_1024_{i}", "status": "STARTING"}
        for i in range(PENDING_DOWNLOADER_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_DOWNLOADER_1024_{i}", "status": "RUNNING"}
        for i in range(RUNNING_DOWNLOADER_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_SALMON_1024_{i}", "status": "SUBMITTED"}
        for i in range(PENDING_SALMON_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_SALMON_1024_{i}", "status": "RUNNING"}
        for i in range(RUNNING_SALMON_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_AFFY_TO_PCL_1024_{i}", "status": "PENDING"}
        for i in range(PENDING_AFFY_JOBS)
    ],
    *[
        {"jobName": f"tests_dev_AFFY_TO_PCL_1024_{i}", "status": "RUNNING"}
        for i in range(RUNNING_AFFY_JOBS)
    ],
]

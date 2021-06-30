##
# Contains DownloaderJobListView, DownloaderJobDetailView and the needed serializer
##

from django.utils.decorators import method_decorator
from rest_framework import filters, generics, serializers

import boto3
from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import DownloaderJob
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)

AWS_REGION = get_env_variable(
    "AWS_REGION", "us-east-1"
)  # Default to us-east-1 if the region variable can't be found

# Job definitons are AWS objects so they have to be namespaced for our stack.
JOB_DEFINITION_PREFIX = get_env_variable("JOB_DEFINITION_PREFIX", "")

batch = boto3.client("batch", region_name=AWS_REGION)


class DownloaderJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = DownloaderJob
        fields = (
            "id",
            "downloader_task",
            "ram_amount",
            "num_retries",
            "retried",
            "was_recreated",
            "worker_id",
            "worker_version",
            "batch_job_id",
            "batch_job_queue",
            "failure_reason",
            "success",
            "original_files",
            "start_time",
            "end_time",
            "created_at",
            "last_modified",
        )
        read_only_fields = fields


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="sample_accession_code",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="List the downloader jobs associated with a sample",
            ),
        ]
    ),
)
class DownloaderJobListView(generics.ListAPIView):
    """
    List of all DownloaderJob
    """

    model = DownloaderJob
    serializer_class = DownloaderJobSerializer
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = DownloaderJobSerializer.Meta.fields
    ordering_fields = ("id", "created_at")
    ordering = ("-id",)

    def list(self, request, *args, **kwargs):
        response = super(DownloaderJobListView, self).list(request, args, kwargs)

        results = response.data["results"]
        batch_job_ids = [job["batch_job_id"] for job in results if job.get("batch_job_id", None)]
        running_job_ids = set()
        if batch_job_ids:
            try:
                described_jobs = batch.describe_jobs(jobs=batch_job_ids)
                for job in described_jobs["jobs"]:
                    if job["status"] in ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]:
                        running_job_ids.add(job["jobId"])
            except Exception as e:
                logger.exception(f"Failure to query about batch_job_ids.")

        for result in results:
            batch_job_id = result.get("batch_job_id", None)
            result["is_queued"] = bool(batch_job_id and batch_job_id in running_job_ids)

        return response

    def get_queryset(self):
        invalid_filters = check_filters(self, ["sample_accession_code"])

        if invalid_filters:
            raise InvalidFilters(invalid_filters=invalid_filters)

        queryset = DownloaderJob.objects.all()

        sample_accession_code = self.request.query_params.get("sample_accession_code", None)
        if sample_accession_code:
            queryset = queryset.filter(
                original_files__samples__accession_code=sample_accession_code
            ).distinct()

        return queryset


class DownloaderJobDetailView(generics.RetrieveAPIView):
    """ Retrieves a DownloaderJob by ID """

    lookup_field = "id"
    model = DownloaderJob
    queryset = DownloaderJob.objects.all()
    serializer_class = DownloaderJobSerializer

    def get(self, request, *args, **kwargs):
        response = super(DownloaderJobDetailView, self).get(request, args, kwargs)

        if "batch_job_id" in response.data and response.data["batch_job_id"]:
            try:
                described_jobs = batch.describe_jobs(jobs=[response.data["batch_job_id"]])
                response.data["is_queued"] = described_jobs["jobs"][0]["status"] in [
                    "SUBMITTED",
                    "PENDING",
                    "RUNNABLE",
                    "STARTING",
                    "RUNNING",
                ]
            except Exception as e:
                logger.exception(
                    f"Failure to query about batch_job_id.",
                    batch_job_id=response.data["batch_job_id"],
                )
        else:
            response.data["is_queued"] = False

        return response

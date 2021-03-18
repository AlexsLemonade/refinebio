##
# Contains DownloaderJobListView, DownloaderJobDetailView and the needed serializer
##

from django.utils.decorators import method_decorator
from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_common.models import DownloaderJob
from data_refinery_common.utils import get_nomad_jobs


class DownloaderJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = DownloaderJob
        fields = (
            "id",
            "downloader_task",
            "num_retries",
            "retried",
            "was_recreated",
            "worker_id",
            "worker_version",
            "batch_job_id",
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
            openapi.Parameter(
                name="nomad",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Only return jobs that are in the nomad queue currently",
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

    def get_queryset(self):
        invalid_filters = check_filters(self, ["sample_accession_code", "nomad"])

        if invalid_filters:
            raise InvalidFilters(invalid_filters=invalid_filters)

        queryset = DownloaderJob.objects.all()

        sample_accession_code = self.request.query_params.get("sample_accession_code", None)
        if sample_accession_code:
            queryset = queryset.filter(
                original_files__samples__accession_code=sample_accession_code
            ).distinct()

        nomad = self.request.query_params.get("nomad", None)
        if nomad:
            running_nomad_jobs_ids = [
                job["ID"] for job in get_nomad_jobs() if job["Status"] == "running"
            ]
            queryset = queryset.filter(batch_job_id__in=running_nomad_jobs_ids)

        return queryset


class DownloaderJobDetailView(generics.RetrieveAPIView):
    """ Retrieves a DownloaderJob by ID """

    lookup_field = "id"
    model = DownloaderJob
    queryset = DownloaderJob.objects.all()
    serializer_class = DownloaderJobSerializer

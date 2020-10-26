##
# Contains OriginalFileListView, OriginalFileDetailView, and needed serializers
##

from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_api.views.relation_serializers import (
    DetailedExperimentSampleSerializer,
    DownloaderJobRelationSerializer,
    ProcessorJobRelationSerializer,
)
from data_refinery_common.models import OriginalFile


class OriginalFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "sha1",
            "samples",
            "processor_jobs",
            "downloader_jobs",
            "source_url",
            "source_filename",
            "is_downloaded",
            "is_archive",
            "has_raw",
            "created_at",
            "last_modified",
        )


class DetailedOriginalFileSerializer(OriginalFileSerializer):
    samples = DetailedExperimentSampleSerializer(many=True)
    processor_jobs = ProcessorJobRelationSerializer(many=True)
    downloader_jobs = DownloaderJobRelationSerializer(many=True)


class OriginalFileListSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "samples",
            "processor_jobs",
            "downloader_jobs",
            "size_in_bytes",
            "sha1",
            "source_url",
            "is_archive",
            "source_filename",
            "has_raw",
            "created_at",
            "last_modified",
        )


class OriginalFileListView(generics.ListAPIView):
    """
    original_files_list

    List Original Files that are associated with Samples. These are the files we proccess.

    """

    queryset = OriginalFile.objects.all()
    serializer_class = OriginalFileListSerializer
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = OriginalFileListSerializer.Meta.fields
    ordering_fields = (
        "id",
        "created_at",
        "last_modified",
    )
    ordering = ("-id",)

    def get_queryset(self):
        invalid_filters = check_filters(self)

        if invalid_filters:
            raise InvalidFilters(invalid_filters)

        return self.queryset


class OriginalFileDetailView(generics.RetrieveAPIView):
    """
    Retrieves an Original File by its ID
    """

    lookup_field = "id"
    queryset = OriginalFile.objects.all()
    serializer_class = DetailedOriginalFileSerializer

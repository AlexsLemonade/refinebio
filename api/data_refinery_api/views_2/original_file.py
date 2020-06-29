from rest_framework import filters, serializers, generics

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_common.models import OriginalFile

# This seems to be unused?
class OriginalFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "sha1",
            "source_url",
            "source_filename",
            "is_downloaded",
            "is_archive",
            "has_raw",
            "created_at",
            "last_modified",
        )


class OriginalFileListSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "samples",
            "size_in_bytes",
            "sha1",
            "samples",
            "processor_jobs",
            "downloader_jobs",
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

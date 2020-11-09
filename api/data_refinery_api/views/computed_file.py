##
# Contains ComputedFileListView, ComputedFileDetailView, and needed serializers
##

from django.core.exceptions import ValidationError
from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_api.views.relation_serializers import (
    ComputationalResultNoFilesRelationSerializer,
    ComputationalResultRelationSerializer,
    DetailedExperimentSampleSerializer,
)
from data_refinery_common.models import APIToken, ComputedFile


class ComputedFileListSerializer(serializers.ModelSerializer):
    result = ComputationalResultNoFilesRelationSerializer(many=False)
    samples = DetailedExperimentSampleSerializer(many=True)
    compendia_organism_name = serializers.CharField(
        source="compendia_organism__name", read_only=True
    )

    def __init__(self, *args, **kwargs):
        super(ComputedFileListSerializer, self).__init__(*args, **kwargs)
        if "context" in kwargs:
            # only include the field `download_url` if a valid token is specified
            # the token lookup happens in the view.
            if "token" not in kwargs["context"]:
                self.fields.pop("download_url")

    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "samples",
            "size_in_bytes",
            "is_qn_target",
            "is_smashable",
            "is_qc",
            "is_compendia",
            "quant_sf_only",
            "compendia_version",
            "compendia_organism_name",
            "sha1",
            "s3_bucket",
            "s3_key",
            "s3_url",
            "download_url",
            "created_at",
            "last_modified",
            "result",
        )
        extra_kwargs = {
            "download_url": {
                "help_text": "This will contain an url to download the file. You must send a valid [token](#tag/token) in order to receive this."
            }
        }


class DetailedComputedFileSerializer(serializers.ModelSerializer):
    result = ComputationalResultRelationSerializer(many=False, read_only=False)
    samples = DetailedExperimentSampleSerializer(many=True)
    compendia_organism_name = serializers.CharField(
        source="compendia_organism__name", read_only=True
    )

    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "samples",
            "size_in_bytes",
            "is_qn_target",
            "is_smashable",
            "is_qc",
            "is_compendia",
            "quant_sf_only",
            "compendia_version",
            "compendia_organism_name",
            "sha1",
            "s3_bucket",
            "s3_key",
            "s3_url",
            "download_url",
            "created_at",
            "last_modified",
            "result",
        )


class ComputedFileListView(generics.ListAPIView):
    """
    computed_files_list

    ComputedFiles are representation of files created by data-refinery processes.

    This can also be used to fetch all the compendia files we have generated with:
    ```
    GET /computed_files?is_compendia=True&is_public=True
    ```
    """

    queryset = ComputedFile.objects.all()
    serializer_class = ComputedFileListSerializer
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = (
        "id",
        "samples",
        "is_qn_target",
        "is_smashable",
        "is_qc",
        "is_compendia",
        "quant_sf_only",
        "svd_algorithm",
        "compendia_version",
        "created_at",
        "last_modified",
        "result__id",
    )
    ordering_fields = (
        "id",
        "created_at",
        "last_modified",
        "compendia_version",
    )
    ordering = ("-id",)

    def get_queryset(self):
        invalid_filters = check_filters(self)

        if invalid_filters:
            raise InvalidFilters(invalid_filters=invalid_filters)

        return self.queryset

    def get_serializer_context(self):
        """
        Extra context provided to the serializer class.
        """
        serializer_context = super(ComputedFileListView, self).get_serializer_context()
        token_id = self.request.META.get("HTTP_API_KEY", None)
        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return {**serializer_context, "token": token}
        except (APIToken.DoesNotExist, ValidationError):
            return serializer_context


class ComputedFileDetailView(generics.RetrieveAPIView):
    """
    Retrieves a computed file by its ID
    """

    lookup_field = "id"
    queryset = ComputedFile.objects.all()
    serializer_class = DetailedComputedFileSerializer

##
# Contains ComputationalResultListView, ComputationalResultDetailView, and needed serializers
##

from django.core.exceptions import ValidationError
from rest_framework import generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.views.relation_serializers import (
    ComputationalResultAnnotationRelationSerializer,
    ComputedFileRelationSerializer,
    ComputedFileWithUrlRelationSerializer,
    OrganismIndexRelationSerializer,
    ProcessorRelationSerializer,
)
from data_refinery_common.models import APIToken, ComputationalResult

##
# Serializers
##


class ComputationalResultSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationRelationSerializer(
        many=True, source="computationalresultannotation_set"
    )
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)
    files = ComputedFileRelationSerializer(many=True, source="computedfile_set")

    class Meta:
        model = ComputationalResult
        fields = (
            "id",
            "commands",
            "processor",
            "is_ccdl",
            "annotations",
            "files",
            "organism_index",
            "time_start",
            "time_end",
            "created_at",
            "last_modified",
        )


class ComputationalResultWithUrlSerializer(ComputationalResultSerializer):
    files = ComputedFileWithUrlRelationSerializer(many=True, source="computedfile_set")


class DetailedComputationalResultSerializer(ComputationalResultSerializer):
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)


class DetailedComputationalResultWithUrlSerializer(ComputationalResultWithUrlSerializer):
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)


##
# Views
##


class ComputationalResultListView(generics.ListAPIView):
    """
    computational_result_list

    This lists all `ComputationalResult`. Each one contains meta-information about the output of a computer process. (Ex Salmon).

    This can return valid S3 urls if a valid [token](#tag/token) is sent in the header `HTTP_API_KEY`.
    """

    queryset = ComputationalResult.public_objects.all()
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ["processor__id"]

    def get_serializer_class(self):
        token_id = self.request.META.get("HTTP_API_KEY", None)

        try:
            # Verify that a token with that id exists
            APIToken.objects.get(id=token_id, is_activated=True)
            return ComputationalResultWithUrlSerializer
        except (APIToken.DoesNotExist, ValidationError):
            return ComputationalResultSerializer

    def filter_queryset(self, queryset):
        filter_dict = self.request.query_params.dict()
        filter_dict.pop("limit", None)
        filter_dict.pop("offset", None)
        return queryset.filter(**filter_dict)


class ComputationalResultDetailView(generics.RetrieveAPIView):
    """
    Retrieves a computational result by its ID
    """

    lookup_field = "id"
    queryset = ComputationalResult.public_objects.all()

    def get_serializer_class(self):
        token_id = self.request.META.get("HTTP_API_KEY", None)

        try:
            # Verify that a token with that id exists
            APIToken.objects.get(id=token_id, is_activated=True)
            return DetailedComputationalResultWithUrlSerializer
        except (APIToken.DoesNotExist, ValidationError):
            return DetailedComputationalResultSerializer

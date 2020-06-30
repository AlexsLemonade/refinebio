##
# Contains List and Detail views for compendium results along with needed serializers
##

from django.db.models.expressions import F, Q
from django.db.models import OuterRef, Subquery
from django.utils.decorators import method_decorator
from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_api.views.organism import OrganismSerializer
from data_refinery_api.views.relation_serializers import (
    ComputationalResultNoFilesRelationSerializer,
)

from data_refinery_common.models import APIToken, CompendiumResult, ComputedFile

from data_refinery_api.views.relation_serializers import (
    ComputedFileRelationSerializer,
    ComputedFileWithUrlRelationSerializer,
)


class CompendiumResultSerializer(serializers.ModelSerializer):
    primary_organism_name = serializers.StringRelatedField(
        read_only=True, source="primary_organism"
    )
    organism_names = serializers.StringRelatedField(many=True, source="organisms", read_only=True)
    computed_file = ComputedFileRelationSerializer(source="get_computed_file", read_only=True)

    class Meta:
        model = CompendiumResult
        fields = (
            "id",
            "primary_organism_name",
            "organism_names",
            "svd_algorithm",
            "quant_sf_only",
            "compendium_version",
            "computed_file",
        )
        read_only_fields = fields


class CompendiumResultWithUrlSerializer(serializers.ModelSerializer):
    primary_organism_name = serializers.StringRelatedField(
        read_only=True, source="primary_organism"
    )
    organism_names = serializers.StringRelatedField(many=True, source="organisms", read_only=True)
    computed_file = ComputedFileWithUrlRelationSerializer(
        source="get_computed_file", read_only=True
    )

    class Meta:
        model = CompendiumResult
        fields = (
            "id",
            "primary_organism_name",
            "organism_names",
            "svd_algorithm",
            "quant_sf_only",
            "compendium_version",
            "computed_file",
        )
        read_only_fields = fields


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="latest_version",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_BOOLEAN,
                description="`True` will only return the highest `compendium_version` for each primary_organism.",
            ),
            openapi.Parameter(
                name="quant_sf_only",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_BOOLEAN,
                description="`True` for RNA-seq Sample Compendium results or `False` for quantile normalized.",
            ),
        ]
    ),
)
class CompendiumResultListView(generics.ListAPIView):
    """
    List all CompendiaResults with filtering.
    """

    model = CompendiumResult
    queryset = CompendiumResult.objects.all()
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = (
        "primary_organism__name",
        "compendium_version",
        "quant_sf_only",
        "result__id",
    )
    ordering_fields = ("primary_organism__name", "compendium_version", "id")
    ordering = ("primary_organism__name",)

    def get_queryset(self):
        public_result_queryset = CompendiumResult.objects.filter(result__is_public=True)
        latest_version = self.request.query_params.get("latest_version", False)
        if latest_version:
            version_filter = Q(
                primary_organism=OuterRef("primary_organism"),
                quant_sf_only=OuterRef("quant_sf_only"),
            )
            latest_version = (
                public_result_queryset.filter(version_filter)
                .order_by("-compendium_version")
                .values("compendium_version")
            )
            return public_result_queryset.annotate(
                latest_version=Subquery(latest_version[:1])
            ).filter(compendium_version=F("latest_version"))

        return public_result_queryset

    def get_serializer_class(self):
        try:
            token_id = self.request.META.get("HTTP_API_KEY", None)
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return CompendiumResultWithUrlSerializer
        except Exception:  # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return CompendiumResultSerializer


class CompendiumResultDetailView(generics.RetrieveAPIView):
    """
    Get a specific Compendium Result
    """

    model = CompendiumResult
    queryset = CompendiumResult.objects.filter(is_public=True)
    lookup_field = "id"

    def get_serializer_class(self):
        try:
            token_id = self.request.META.get("HTTP_API_KEY", None)
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return CompendiumResultWithUrlSerializer
        except Exception:  # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return CompendiumResultSerializer

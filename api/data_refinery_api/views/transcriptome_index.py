##
# Contains TranscriptomeIndexListView, TranscriptomeIndexDetailView, and needed serializer
##

from django.utils.decorators import method_decorator
from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_common.models import OrganismIndex


class OrganismIndexSerializer(serializers.ModelSerializer):

    organism_name = serializers.StringRelatedField(source="organism", read_only=True)
    download_url = serializers.SerializerMethodField()

    class Meta:
        model = OrganismIndex

        fields = (
            "id",
            "assembly_name",
            "organism_name",
            "source_version",
            "index_type",
            "salmon_version",
            "download_url",
            "result_id",
            "last_modified",
        )
        read_only_fields = fields

    def get_download_url(self, obj):
        computed_file = obj.get_computed_file()
        if computed_file is not None:
            return computed_file.s3_url
        return None


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="organism_name",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Organism name. Eg. `MUS_MUSCULUS`",
            ),
            openapi.Parameter(
                name="length",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Short hand for `index_type` Eg. `short` or `long`",
            ),
            openapi.Parameter(
                name="salmon_version",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Eg. `salmon 0.13.1`",
            ),
            openapi.Parameter(
                name="index_type",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Eg. `TRANSCRIPTOME_LONG`",
            ),
        ]
    ),
)
class TranscriptomeIndexListView(generics.ListAPIView):
    """
    List all Transcriptome Indices. These are a special type of process result,
    necessary for processing other SRA samples.
    """

    serializer_class = OrganismIndexSerializer
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = ["salmon_version", "index_type"]
    ordering_fields = ("created_at", "salmon_version")
    ordering = ("-created_at",)

    def get_queryset(self):
        queryset = OrganismIndex.public_objects.all()

        organism_name = self.request.query_params.get("organism_name", None)
        if organism_name is not None:
            queryset = queryset.filter(organism__name=organism_name.upper())

        # https://github.com/AlexsLemonade/refinebio/issues/2459
        # It looks like when we set `result_id` as a filterset field,
        # django_forms goes nuts and tries to call __str__ on every single
        # computational result in our database trying to find all of the
        # different possible computational_results. So let's just take care of
        # this one ourselves.
        result_id = self.request.query_params.get("result_id", None)
        if result_id is not None:
            queryset = queryset.filter(result_id=result_id)

        length = self.request.query_params.get("length", None)
        if length is not None:
            index_type = "TRANSCRIPTOME_{}".format(length.upper())
            queryset = queryset.filter(index_type=index_type)

        return queryset


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="id",
                in_=openapi.IN_PATH,
                type=openapi.TYPE_NUMBER,
                description="Transcriptome Index Id eg `1`",
            ),
        ]
    ),
)
class TranscriptomeIndexDetailView(generics.RetrieveAPIView):
    """
    Gets the S3 url associated with the organism and length, along with other metadata about
    the transcriptome index we have stored.
    """

    serializer_class = OrganismIndexSerializer
    lookup_field = "id"
    queryset = OrganismIndex.public_objects.all()

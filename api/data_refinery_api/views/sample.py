##
# Contains SampleListView, SampleDetailView, and needed serializers
##

from django.db.models import Prefetch
from django.db.models.expressions import Q
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_api.views.relation_serializers import (
    OrganismIndexRelationSerializer,
    OrganismRelationSerializer,
    ProcessorRelationSerializer,
)
from data_refinery_common.models import (
    ComputationalResult,
    Dataset,
    Experiment,
    Sample,
    SampleAnnotation,
)


class SampleAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleAnnotation
        fields = (
            "data",
            "is_ccdl",
            "created_at",
            "last_modified",
        )


class DetailedSamplesComputationalResultSerializer(serializers.ModelSerializer):
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)

    class Meta:
        model = ComputationalResult
        fields = (
            "id",
            "processor",
            "organism_index",
        )


class DetailedSampleSerializer(serializers.ModelSerializer):
    annotations = SampleAnnotationSerializer(many=True, source="sampleannotation_set")
    organism = OrganismRelationSerializer(many=False)
    results = DetailedSamplesComputationalResultSerializer(many=True)

    class Meta:
        model = Sample
        fields = (
            "id",
            "title",
            "accession_code",
            "source_database",
            "organism",
            "platform_accession_code",
            "platform_name",
            "pretty_platform",
            "technology",
            "manufacturer",
            "protocol_info",
            "annotations",
            "results",
            "source_archive_url",
            "has_raw",
            "sex",
            "age",
            "specimen_part",
            "genotype",
            "disease",
            "disease_stage",
            "cell_line",
            "treatment",
            "race",
            "subject",
            "compound",
            "time",
            "is_processed",
            "created_at",
            "last_modified",
            "original_files",
            "computed_files",
            "experiment_accession_codes",
        )


@method_decorator(
    name="get",
    decorator=swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="dataset_id",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Filters the result and only returns samples that are added to a dataset.",
            ),
            openapi.Parameter(
                name="experiment_accession_code",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Filters the result and only returns only the samples associated with an experiment accession code.",
            ),
            openapi.Parameter(
                name="accession_codes",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Provide a list of sample accession codes separated by commas and the endpoint will only return information about these samples.",
            ),
        ]
    ),
)
class SampleListView(generics.ListAPIView):
    """ Returns detailed information about Samples """

    model = Sample
    serializer_class = DetailedSampleSerializer
    filter_backends = (filters.OrderingFilter, DjangoFilterBackend)
    ordering_fields = "__all__"
    ordering = "-is_processed"
    filterset_fields = (
        "title",
        "organism",
        "source_database",
        "source_archive_url",
        "has_raw",
        "platform_name",
        "technology",
        "manufacturer",
        "sex",
        "age",
        "specimen_part",
        "genotype",
        "disease",
        "disease_stage",
        "cell_line",
        "treatment",
        "race",
        "subject",
        "compound",
        "time",
        "is_processed",
        "is_public",
    )

    def get_queryset(self):
        """
        ref https://www.django-rest-framework.org/api-guide/filtering/#filtering-against-query-parameters
        """
        invalid_filters = check_filters(
            self,
            special_filters=[
                "ids",
                "organism__name",
                "dataset_id",
                "experiment_accession_code",
                "accession_codes",
            ],
        )

        if invalid_filters:
            raise InvalidFilters(invalid_filters)

        queryset = (
            Sample.public_objects.prefetch_related("organism")
            .prefetch_related(
                Prefetch("results", queryset=ComputationalResult.objects.order_by("time_start"))
            )
            .prefetch_related("results__processor")
            .prefetch_related("results__computationalresultannotation_set")
            .prefetch_related("results__computedfile_set")
            .filter(**self.get_query_params_filters())
        )

        # case insensitive search https://docs.djangoproject.com/en/2.1/ref/models/querysets/#icontains
        filter_by = self.request.query_params.get("filter_by", None)
        if filter_by:
            queryset = queryset.filter(
                Q(accession_code__icontains=filter_by)
                | Q(title__icontains=filter_by)
                | Q(sex__icontains=filter_by)
                | Q(age__icontains=filter_by)
                | Q(specimen_part__icontains=filter_by)
                | Q(genotype__icontains=filter_by)
                | Q(disease__icontains=filter_by)
                | Q(disease_stage__icontains=filter_by)
                | Q(cell_line__icontains=filter_by)
                | Q(treatment__icontains=filter_by)
                | Q(race__icontains=filter_by)
                | Q(subject__icontains=filter_by)
                | Q(compound__icontains=filter_by)
                | Q(time__icontains=filter_by)
            )

        return queryset

    def get_query_params_filters(self):
        """ We do advanced filtering on the queryset depending on the query parameters.
            This returns the parameters that should be used for that. """
        filter_dict = dict()

        ids = self.request.query_params.get("ids", None)
        if ids is not None:
            ids = [int(x) for x in ids.split(",")]
            filter_dict["pk__in"] = ids

        experiment_accession_code = self.request.query_params.get("experiment_accession_code", None)
        if experiment_accession_code:
            experiment = get_object_or_404(
                Experiment.objects.values("id"), accession_code=experiment_accession_code
            )
            filter_dict["experiments__in"] = [experiment["id"]]

        accession_codes = self.request.query_params.get("accession_codes", None)
        if accession_codes:
            accession_codes = accession_codes.split(",")
            filter_dict["accession_code__in"] = accession_codes

        dataset_id = self.request.query_params.get("dataset_id", None)
        if dataset_id:
            dataset = get_object_or_404(Dataset, id=dataset_id)
            # Python doesn't provide a prettier way of doing this that I know about.
            filter_dict["accession_code__in"] = [
                item for sublist in dataset.data.values() for item in sublist
            ]

        # Accept Organism in both name and ID form
        organism_name = self.request.query_params.get("organism__name", None)
        if organism_name:
            filter_dict["organism__name"] = organism_name

        return filter_dict


class SampleDetailView(generics.RetrieveAPIView):
    """ Retrieve the details for a Sample given its accession code """

    lookup_field = "accession_code"
    queryset = Sample.public_objects.all()
    serializer_class = DetailedSampleSerializer

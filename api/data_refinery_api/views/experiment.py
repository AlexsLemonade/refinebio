##
# Contains ExperimentListView, ExperimentDetailView, and needed serializers
##

from django.db.models import Count, Q
from rest_framework import generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_api.views.relation_serializers import DetailedExperimentSampleSerializer
from data_refinery_common.models import Experiment, ExperimentAnnotation


class ExperimentAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExperimentAnnotation
        fields = (
            "data",
            "is_ccdl",
            "created_at",
            "last_modified",
        )


class ExperimentSerializer(serializers.ModelSerializer):
    organism_names = serializers.StringRelatedField(many=True, source="organisms", read_only=True)
    platforms = serializers.ReadOnlyField()
    processed_samples = serializers.StringRelatedField(many=True)
    total_samples_count = serializers.IntegerField(read_only=True)
    sample_metadata = serializers.ReadOnlyField(source="get_sample_metadata_fields")
    technologies = serializers.ReadOnlyField(source="get_sample_technologies")
    pretty_platforms = serializers.ReadOnlyField()

    class Meta:
        model = Experiment
        fields = (
            "id",
            "title",
            "description",
            "accession_code",
            "alternate_accession_code",
            "source_database",
            "source_url",
            "platforms",
            "pretty_platforms",
            "processed_samples",
            "has_publication",
            "publication_title",
            "publication_doi",
            "publication_authors",
            "pubmed_id",
            "total_samples_count",
            "organism_names",
            "submitter_institution",
            "created_at",
            "last_modified",
            "source_first_published",
            "source_last_modified",
            "sample_metadata",
            "technologies",
        )

    @staticmethod
    def setup_eager_loading(queryset):
        """ Perform necessary eager loading of data. """
        queryset = queryset.prefetch_related("samples").prefetch_related("organisms")

        # Multiple count annotations
        queryset = queryset.annotate(
            total_samples_count=Count("samples", unique=True),
            processed_samples_count=Count("samples", filter=Q(samples__is_processed=True)),
        )

        return queryset


class DetailedExperimentSerializer(serializers.ModelSerializer):
    annotations = ExperimentAnnotationSerializer(many=True, source="experimentannotation_set")
    samples = DetailedExperimentSampleSerializer(many=True)
    sample_metadata = serializers.ReadOnlyField(source="sample_metadata_fields")
    organism_names = serializers.StringRelatedField(many=True, source="organisms", read_only=True)

    class Meta:
        model = Experiment
        fields = (
            "id",
            "title",
            "description",
            "annotations",
            "samples",
            "protocol_description",
            "accession_code",
            "alternate_accession_code",
            "source_database",
            "source_url",
            "has_publication",
            "publication_title",
            "publication_doi",
            "publication_authors",
            "pubmed_id",
            "source_first_published",
            "source_last_modified",
            "submitter_institution",
            "last_modified",
            "created_at",
            "organism_names",
            "sample_metadata",
            "num_total_samples",
            "num_processed_samples",
            "num_downloadable_samples",
        )


class ExperimentListView(generics.ListAPIView):
    """ Paginated list of all experiments. Advanced filtering can be done with the `/search` endpoint. """

    model = Experiment
    queryset = Experiment.public_objects.all()
    serializer_class = ExperimentSerializer
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = (
        "title",
        "description",
        "accession_code",
        "alternate_accession_code",
        "source_database",
        "source_url",
        "has_publication",
        "publication_title",
        "publication_doi",
        "pubmed_id",
        "organisms",
        "submitter_institution",
        "created_at",
        "last_modified",
        "source_first_published",
        "source_last_modified",
    )

    def get_queryset(self):
        invalid_filters = check_filters(self)

        if invalid_filters:
            raise InvalidFilters(invalid_filters)

        return self.queryset


class ExperimentDetailView(generics.RetrieveAPIView):
    """ Retrieve details for an experiment given it's accession code """

    lookup_field = "accession_code"
    queryset = Experiment.public_objects.all()
    serializer_class = DetailedExperimentSerializer

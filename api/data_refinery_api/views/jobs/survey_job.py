##
# Contains SurveyJobListView, SurveyJobDetailView, and the needed serializer
##

from rest_framework import filters, generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_common.models import SurveyJob


class SurveyJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = SurveyJob
        fields = (
            "id",
            "source_type",
            "success",
            "start_time",
            "end_time",
            "created_at",
            "last_modified",
        )


class SurveyJobListView(generics.ListAPIView):
    """
    List of all SurveyJob.
    """

    model = SurveyJob
    queryset = SurveyJob.objects.all()
    serializer_class = SurveyJobSerializer
    filter_backends = (
        DjangoFilterBackend,
        filters.OrderingFilter,
    )
    filterset_fields = SurveyJobSerializer.Meta.fields
    ordering_fields = ("id", "created_at")
    ordering = ("-id",)

    def get_queryset(self):
        invalid_filters = check_filters(self)

        if invalid_filters:
            raise InvalidFilters(invalid_filters=invalid_filters)

        return self.queryset


class SurveyJobDetailView(generics.RetrieveAPIView):
    """ Retrieves a SurveyJob by ID """

    lookup_field = "id"
    model = SurveyJob
    queryset = SurveyJob.objects.all()
    serializer_class = SurveyJobSerializer

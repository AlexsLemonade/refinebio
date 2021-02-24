from rest_framework import generics, serializers

from django_filters.rest_framework import DjangoFilterBackend

from data_refinery_api.exceptions import InvalidFilters
from data_refinery_api.utils import check_filters
from data_refinery_common.models import Organism


class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = (
            "name",
            "taxonomy_id",
            "has_compendia",
            "has_quantfile_compendia",
        )
        read_only_fields = fields


class OrganismListView(generics.ListAPIView):
    """
    Paginated list of all the available organisms.
    """

    filter_backends = (DjangoFilterBackend,)
    filterset_fields = ("has_compendia", "has_quantfile_compendia")
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer

    def get_queryset(self):
        invalid_filters = check_filters(self)

        if invalid_filters:
            raise InvalidFilters(invalid_filters=invalid_filters)

        return self.queryset


class OrganismDetailView(generics.RetrieveAPIView):
    """
    Retrieves an organism by its name.
    """

    lookup_field = "name"
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer

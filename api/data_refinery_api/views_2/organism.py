from rest_framework import serializers, generics
from data_refinery_common.models import Organism


class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = (
            "name",
            "taxonomy_id",
        )


class OrganismListView(generics.ListAPIView):
    """
    Paginated list of all the available organisms.
    """

    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer

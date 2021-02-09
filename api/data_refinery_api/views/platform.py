from rest_framework import generics, serializers

from data_refinery_common.models import Sample


class PlatformSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sample
        fields = (
            "platform_accession_code",
            "platform_name",
        )
        read_only_fields = fields


class PlatformListView(generics.ListAPIView):
    """
    Unpaginated list of all the available "platform" information
    """

    serializer_class = PlatformSerializer
    paginator = None

    def get_queryset(self):
        return (
            Sample.public_objects.all()
            .values("platform_accession_code", "platform_name")
            .distinct()
        )

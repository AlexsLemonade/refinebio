from rest_framework import generics, serializers

from data_refinery_common.models import Processor


class ProcessorSerializer(serializers.ModelSerializer):
    class Meta:
        model = Processor
        fields = ("id", "name", "version", "docker_image", "environment")
        read_only_fields = fields


class ProcessorListView(generics.ListAPIView):
    """List all processors."""

    queryset = Processor.objects.all()
    serializer_class = ProcessorSerializer


class ProcessorDetailView(generics.RetrieveAPIView):
    """Retrieves a processor by its ID"""

    lookup_field = "id"
    queryset = Processor.objects.all()
    serializer_class = ProcessorSerializer

from rest_framework import serializers, generics
from data_refinery_common.models import Processor


class ProcessorSerializer(serializers.ModelSerializer):
    class Meta:
        model = Processor
        fields = ("id", "name", "version", "docker_image", "environment")


class ProcessorListView(generics.ListAPIView):
    """List all processors."""

    queryset = Processor.objects.all()
    serializer_class = ProcessorSerializer

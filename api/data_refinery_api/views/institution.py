from rest_framework import generics, serializers

from data_refinery_common.models import Experiment


class InstitutionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ("submitter_institution",)
        read_only_fields = fields


class InstitutionListView(generics.ListAPIView):
    """
    Unpaginated list of all the available "institution" information
    """

    serializer_class = InstitutionSerializer
    paginator = None

    def get_queryset(self):
        return Experiment.public_objects.all().values("submitter_institution").distinct()

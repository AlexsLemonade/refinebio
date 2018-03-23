from django.http import Http404
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from data_refinery_common.models import Experiment
from data_refinery_api.serializers import ExperimentSerializer

class ExperimentList(APIView):
    """
    List all Experiments
    """
    def get(self, request, format=None):
        experiments = Experiment.objects.all()
        serializer = ExperimentSerializer(experiments, many=True)
        return Response(serializer.data)

class ExperimentDetail(APIView):
    """
    Retriev an Experiment instance.
    """
    def get_object(self, pk):
        try:
            return Experiment.objects.get(pk=pk)
        except Experiment.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        experiment = self.get_object(pk)
        serializer = ExperimentSerializer(experiment)
        return Response(serializer.data)
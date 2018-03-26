from django.conf import settings
from rest_framework.settings import api_settings
from django.http import Http404

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from data_refinery_common.models import Experiment, Sample, Organism
from data_refinery_api.serializers import ( 
	ExperimentSerializer, 
	DetailedExperimentSerializer,
	SampleSerializer, 
	DetailedSampleSerializer,
	OrganismSerializer,
	PlatformSerializer,
	InstitutionSerializer
)

class PaginatedAPIView(APIView):
    pagination_class = api_settings.DEFAULT_PAGINATION_CLASS

    @property
    def paginator(self):
        """
        The paginator instance associated with the view, or `None`.
        """
        if not hasattr(self, '_paginator'):
            if self.pagination_class is None:
                self._paginator = None
            else:
                self._paginator = self.pagination_class()
        return self._paginator

    def paginate_queryset(self, queryset):
        """
        Return a single page of results, or `None` if pagination is disabled.
        """
        if self.paginator is None:
            return None
        return self.paginator.paginate_queryset(queryset, self.request, view=self)

    def get_paginated_response(self, data):
        """
        Return a paginated style `Response` object for the given output data.
        """
        assert self.paginator is not None
        return self.paginator.get_paginated_response(data)

##
# Experiments
##

class ExperimentList(PaginatedAPIView):
    """
    List all Experiments
    """

    def get(self, request, format=None):
        experiments = Experiment.objects.all()

        page = self.paginate_queryset(experiments)
        if page is not None:
            serializer = ExperimentSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
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
        serializer = DetailedExperimentSerializer(experiment)
        return Response(serializer.data)

##
# Experiments
##

class SampleList(PaginatedAPIView):
    """
    List all Samples
    """

    def get(self, request, format=None):
        samples = Sample.objects.all()

        page = self.paginate_queryset(samples)
        if page is not None:
            serializer = SampleSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = SampleSerializer(samples, many=True)
            return Response(serializer.data)

class SampleDetail(APIView):
    """
    Retriev an Sample instance.
    """
    def get_object(self, pk):
        try:
            return Sample.objects.get(pk=pk)
        except Sample.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sample = self.get_object(pk)
        serializer = DetailedSampleSerializer(sample)
        return Response(serializer.data)

##
# Search Filter Models
##

class OrganismList(APIView):
    """
	Unpaginated list of all the available organisms
	"""
     
    def get(self, request, format=None):
        organisms = Organism.objects.all()
        serializer = OrganismSerializer(organisms, many=True)
        return Response(serializer.data)

class PlatformList(APIView):
    """
	Unpaginated list of all the available "platform" information
	"""
     
    def get(self, request, format=None):
        experiments = Experiment.objects.all().values("platform_accession_code", "platform_name").distinct()
        serializer = PlatformSerializer(experiments, many=True)
        return Response(serializer.data)

class InstitutionList(APIView):
    """
	Unpaginated list of all the available "institution" information
	"""
     
    def get(self, request, format=None):
        experiments = Experiment.objects.all().values("submitter_institution").distinct()
        serializer = InstitutionSerializer(experiments, many=True)
        return Response(serializer.data)

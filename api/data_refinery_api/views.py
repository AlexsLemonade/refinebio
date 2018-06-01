from django.conf import settings
from django.db.models import Count
from django.db.models.aggregates import Avg
from django.db.models.expressions import F
from django.http import Http404

from django_filters.rest_framework import DjangoFilterBackend

from rest_framework.pagination import LimitOffsetPagination
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.settings import api_settings
from rest_framework import status, filters, generics

from data_refinery_common.models import (
    Experiment,
    Sample,
    Organism,
    ComputationalResult,
    DownloaderJob,
    SurveyJob,
    ProcessorJob,
    Dataset
)
from data_refinery_api.serializers import (
    ExperimentSerializer,
    DetailedExperimentSerializer,
    SampleSerializer,
    DetailedSampleSerializer,
    OrganismSerializer,
    PlatformSerializer,
    InstitutionSerializer,
    ComputationalResultSerializer,

    # Jobs
    SurveyJobSerializer,
    DownloaderJobSerializer,
    ProcessorJobSerializer,

    # Dataset
    CreateDatasetSerializer,
    DatasetSerializer
)

##
# Custom Views
##

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
# Search and Filter
##

# ListAPIView is read-only!
class SearchAndFilter(generics.ListAPIView):
    """
    Search and filter for experiments and samples.

    Ex: search/?search=human&has_publication=True

    """

    queryset = Experiment.objects.all()
    serializer_class = ExperimentSerializer
    pagination_class = LimitOffsetPagination

    filter_backends = (DjangoFilterBackend, filters.SearchFilter,)
    search_fields = ('title', '@description')
    filter_fields = ('has_publication', 'submitter_institution', 'technology', 'source_first_published')


    def list(self, request, *args, **kwargs):
        """ Adds counts on certain filter fields to result JSON."""
        response = super(SearchAndFilter, self).list(request, args, kwargs)

        response.data['filters'] = {}
        response.data['filters']['technology'] = {}
        response.data['filters']['publication'] = {}
        response.data['filters']['organism'] = {}

        qs = self.filter_queryset(self.get_queryset())

        techs = qs.values('technology').annotate(Count('technology', unique=True))
        for tech in techs:
            response.data['filters']['technology'][tech['technology']] = tech['technology__count']

        pubs = qs.values('has_publication').annotate(Count('has_publication', unique=True))
        for pub in pubs:
            if pub['has_publication']:
                response.data['filters']['publication']['has_publication'] = pub['has_publication__count']
        if 'has_publication' not in response.data['filters']['publication']:
            response.data['filters']['publication']['has_publication'] = 0

        organisms = qs.values('organisms__name').annotate(Count('organisms__name', unique=True))
        for organism in organisms:
            response.data['filters']['organism'][organism['organisms__name']] = organism['organisms__name__count']

        return response

##
# Dataset
##

class CreateDatasetView(generics.CreateAPIView):
    """ Creates and returns new Dataset. """

    queryset = Dataset.objects.all()
    serializer_class = CreateDatasetSerializer

class DatasetView(generics.RetrieveUpdateAPIView):
    """ View and modify a single Dataset. Set `start` to `true` to begin smashing and delivery."""

    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer
    lookup_field = 'id'

    def perform_update(self, serializer):
        """ If `start` is set, fire off the job. Disables dataset data updates after that. """
        old_object = self.get_object()
        old_data = old_object.data
        old_aggregate = old_object.aggregate_by
        already_processing = old_object.is_processing
        new_data = serializer.validated_data

        if new_data.get('start'):
            if not already_processing:
                # TODO: Fire off the Smasher job here
                serializer.validated_data['is_processing'] = True
                obj = serializer.save()
                return obj

        # Don't allow critical data updates to jobs that have already been submitted,
        # but do allow email address updating.
        if already_processing:
            serializer.validated_data['data'] = old_data
            serializer.validated_data['aggregate_by'] = old_aggregate
        serializer.save()

##
# Experiments
##

class ExperimentList(PaginatedAPIView):
    """
    List all Experiments.

    Append the pk to the end of this URL to see a detail view.
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        experiments = Experiment.objects.filter(**filter_dict)

        page = self.paginate_queryset(experiments)
        if page is not None:
            serializer = ExperimentSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = ExperimentSerializer(experiments, many=True)
            return Response(serializer.data)

class ExperimentDetail(APIView):
    """
    Retrieve an Experiment instance.
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
# Samples
##

class SampleList(PaginatedAPIView):
    """
    List all Samples.

    Pass in a list of pk to an ids query parameter to filter by id.

    Append the pk to the end of this URL to see a detail view.

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        order_by = filter_dict.pop('order_by', None)
        ids = filter_dict.pop('ids', None)

        if ids is not None:
            ids = [ int(x) for x in ids.split(',')]
            filter_dict['pk__in'] = ids

        samples = Sample.objects.filter(**filter_dict)
        if order_by:
            samples = samples.order_by(order_by)

        page = self.paginate_queryset(samples)
        if page is not None:
            serializer = DetailedSampleSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = DetailedSampleSerializer(samples, many=True)
            return Response(serializer.data)

class SampleDetail(APIView):
    """
    Retrieve a Sample instance.
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
# Results
##

class ResultsList(PaginatedAPIView):
    """
    List all ComputationalResults.

    Append the pk to the end of this URL to see a detail view.

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        results = ComputationalResult.objects.filter(**filter_dict)

        page = self.paginate_queryset(results)
        if page is not None:
            serializer = ComputationalResultSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = ComputationalResultSerializer(results, many=True)
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

##
# Jobs
##

class SurveyJobList(PaginatedAPIView):
    """
    List of all SurveyJob.

	Ex:
	  - ?start_time__lte=2018-03-23T15:29:40.848381Z
	  - ?start_time__lte=2018-03-23T15:29:40.848381Z&start_time__gte=2018-03-23T14:29:40.848381Z
	  - ?success=True

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        jobs = SurveyJob.objects.filter(**filter_dict)

        page = self.paginate_queryset(jobs)
        if page is not None:
            serializer = SurveyJobSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = SurveyJobSerializer(jobs, many=True)
            return Response(serializer.data)

class DownloaderJobList(PaginatedAPIView):
    """
    List of all DownloaderJob
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        jobs = DownloaderJob.objects.filter(**filter_dict)

        page = self.paginate_queryset(jobs)
        if page is not None:
            serializer = DownloaderJobSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = DownloaderJobSerializer(jobs, many=True)
            return Response(serializer.data)

class ProcessorJobList(PaginatedAPIView):
    """
    List of all ProcessorJobs
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        jobs = ProcessorJob.objects.filter(**filter_dict)

        page = self.paginate_queryset(jobs)
        if page is not None:
            serializer = ProcessorJobSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            serializer = ProcessorJobSerializer(jobs, many=True)
            return Response(serializer.data)

###
# Statistics
###

class Stats(APIView):
    """
    Statistics about the health of the system.
    """

    def get(self, request, format=None):
        data = {}
        data['survey_jobs'] = {}
        data['survey_jobs']['total'] = SurveyJob.objects.count()
        data['survey_jobs']['pending'] = SurveyJob.objects.filter(start_time__isnull=True).count()
        data['survey_jobs']['completed'] = SurveyJob.objects.filter(end_time__isnull=False).count()
        data['survey_jobs']['open'] = SurveyJob.objects.filter(start_time__isnull=False, end_time__isnull=True).count()
        # via https://stackoverflow.com/questions/32520655/get-average-of-difference-of-datetime-fields-in-django
        data['survey_jobs']['average_time'] = SurveyJob.objects.filter(start_time__isnull=False, end_time__isnull=False).aggregate(average_time=Avg(F('end_time') - F('start_time')))['average_time']

        data['downloader_jobs'] = {}
        data['downloader_jobs']['total'] = DownloaderJob.objects.count()
        data['downloader_jobs']['pending'] = DownloaderJob.objects.filter(start_time__isnull=True).count()
        data['downloader_jobs']['completed'] = DownloaderJob.objects.filter(end_time__isnull=False).count()
        data['downloader_jobs']['open'] = DownloaderJob.objects.filter(start_time__isnull=False, end_time__isnull=True).count()
        data['downloader_jobs']['average_time'] = DownloaderJob.objects.filter(start_time__isnull=False, end_time__isnull=False).aggregate(average_time=Avg(F('end_time') - F('start_time')))['average_time']

        data['processor_jobs'] = {}
        data['processor_jobs']['total'] = ProcessorJob.objects.count()
        data['processor_jobs']['pending'] = ProcessorJob.objects.filter(start_time__isnull=True).count()
        data['processor_jobs']['completed'] = ProcessorJob.objects.filter(end_time__isnull=False).count()
        data['processor_jobs']['open'] = ProcessorJob.objects.filter(start_time__isnull=False, end_time__isnull=True).count()
        data['processor_jobs']['average_time'] = ProcessorJob.objects.filter(start_time__isnull=False, end_time__isnull=False).aggregate(average_time=Avg(F('end_time') - F('start_time')))['average_time']

        return Response(data)


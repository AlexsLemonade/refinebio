from django.conf import settings
from django.db.models import Count
from django.db.models.aggregates import Avg
from django.db.models.expressions import F
from django.http import Http404, HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404

from django_filters.rest_framework import DjangoFilterBackend
import django_filters

from rest_framework.exceptions import APIException
from rest_framework.pagination import LimitOffsetPagination
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.settings import api_settings
from rest_framework.exceptions import ValidationError
from rest_framework import status, filters, generics

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    Experiment,
    Sample,
    Organism,
    Processor,
    ComputationalResult,
    DownloaderJob,
    SurveyJob,
    ProcessorJob,
    Dataset,
    APIToken,
    ProcessorJobDatasetAssociation,
    OrganismIndex
)
from data_refinery_api.serializers import (
    ExperimentSerializer,
    DetailedExperimentSerializer,
    SampleSerializer,
    DetailedSampleSerializer,
    OrganismSerializer,
    OrganismIndexSerializer,
    PlatformSerializer,
    InstitutionSerializer,
    ComputationalResultSerializer,
    ProcessorSerializer,

    # Job
    SurveyJobSerializer,
    DownloaderJobSerializer,
    ProcessorJobSerializer,

    # Dataset
    CreateDatasetSerializer,
    DatasetSerializer,
    APITokenSerializer
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

class ExperimentFilter(django_filters.FilterSet):
    queryset = Experiment.processed_public_objects.all()
    has_publication = django_filters.BooleanFilter(field_name="has_publication")
    submitter_institution = \
        django_filters.ModelMultipleChoiceFilter(field_name="submitter_institution",
                                                 to_field_name="submitter_institution",
                                                 queryset=queryset)
    submitter_institution.always_filter = False
    technology = django_filters.ModelMultipleChoiceFilter(field_name="technology",
                                                          to_field_name="technology",
                                                          queryset=queryset)
    technology.always_filter = False
    source_first_published = django_filters.DateTimeFilter(field_name="source_first_published")
    organisms__name = django_filters.ModelMultipleChoiceFilter(field_name="organisms__name",
                                                               to_field_name="name",
                                                               queryset=Organism.objects.all())
    organisms__name.always_filter = False
    samples__platform_accession_code = \
        django_filters.ModelMultipleChoiceFilter(field_name="smaples__platform_accession_code",
                                                 to_field_name="platform_accession_code",
                                                 queryset=Sample.objects.all())
    samples__platform_accession_code.always_filter = False

    class Meta:
        model = Experiment
        fields = ['has_publication', 
                        'submitter_institution', 
                        'technology',
                        'source_first_published', 
                        'organisms__name',
                        'samples__platform_accession_code']

# ListAPIView is read-only!
class SearchAndFilter(generics.ListAPIView):
    """
    Search and filter for experiments and samples.

    Ex: search/?search=human&has_publication=True

    """

    queryset = Experiment.processed_public_objects.all()

    serializer_class = ExperimentSerializer
    pagination_class = LimitOffsetPagination

    filter_backends = (DjangoFilterBackend, filters.SearchFilter,)
    filter_class = ExperimentFilter

    # via http://www.django-rest-framework.org/api-guide/filtering/#searchfilter
    # '^' Starts-with search.
    # '=' Exact matches.
    # '@' Full-text search.
    # '$' Regex search.
    search_fields = (   'title',
                        '@description',
                        '@accession_code',
                        '@protocol_description',
                        '@publication_title',
                        'publication_doi',
                        'publication_authors',
                        'pubmed_id',
                        '@submitter_institution',
                        'experimentannotation__data'
                    )

    def list(self, request, *args, **kwargs):
        """ Adds counts on certain filter fields to result JSON."""
        response = super(SearchAndFilter, self).list(request, args, kwargs)

        response.data['filters'] = {}
        response.data['filters']['technology'] = {}
        response.data['filters']['publication'] = {}
        response.data['filters']['organism'] = {}

        qs = self.search_queryset(self.get_queryset())
        techs = qs.values('technology').annotate(Count('technology', unique=True))
        for tech in techs:
            if not tech['technology'] or not tech['technology'].strip():
                continue
            response.data['filters']['technology'][tech['technology']] = tech['technology__count']

        pubs = qs.values('has_publication').annotate(Count('has_publication', unique=True))
        for pub in pubs:
            if pub['has_publication']:
                response.data['filters']['publication']['has_publication'] = pub['has_publication__count']

        organisms = qs.values('organisms__name').annotate(Count('organisms__name', unique=True))
        for organism in organisms:

            # This experiment has no ExperimentOrganism-association, which is bad.
            # This information may still live on the samples though.
            if not organism['organisms__name']:
                continue

            response.data['filters']['organism'][organism['organisms__name']] = organism['organisms__name__count']

        return response

    # We want to determine filters based off of the search term but not the filters to allow for
    # multiple filters of the same type.
    def search_queryset(self, queryset):
        """ Filters the queryset based off of the search term (but not the filters) """
        return filters.SearchFilter().filter_queryset(self.request, queryset, view=self)

##
# Dataset
##

class CreateDatasetView(generics.CreateAPIView):
    """ Creates and returns new Dataset. """

    queryset = Dataset.objects.all()
    serializer_class = CreateDatasetSerializer

class DatasetView(generics.RetrieveUpdateAPIView):
    """ View and modify a single Dataset. Set `start` to `true` along with a valid
    activated API token (from /token/) to begin smashing and delivery.
    """

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

            # Make sure we have a valid activated token.
            token_id = self.request.data.get('token_id')
            try:
                token = APIToken.objects.get(id=token_id, is_activated=True)
            except Exception: # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
                raise APIException("You must provide an active API token ID")

            if not already_processing:

                # Create and dispatch the new job.
                processor_job = ProcessorJob()
                processor_job.pipeline_applied = "SMASHER"
                processor_job.ram_amount = 4096
                processor_job.save()

                pjda = ProcessorJobDatasetAssociation()
                pjda.processor_job = processor_job
                pjda.dataset = old_object
                pjda.save()

                # Hidden method of non-dispatching for testing purposes.
                if not self.request.data.get('no_send_job', False):
                    send_job(ProcessorPipeline.SMASHER, processor_job)

                serializer.validated_data['is_processing'] = True
                obj = serializer.save()
                return obj

        # Don't allow critical data updates to jobs that have already been submitted,
        # but do allow email address updating.
        if already_processing:
            serializer.validated_data['data'] = old_data
            serializer.validated_data['aggregate_by'] = old_aggregate
        serializer.save()

class DatasetStatsView(APIView):
    """ Get stats for a given dataset. Ex:

    {
        "HOMO_SAPIENS": {
            "num_experiments": 5,
            "num_samples": 55 },
        "GALLUS_GALLUS": {
            "num_experiments": 5,
            "num_samples": 55 },
    }

    """

    def get(self, request, id):
        
        dataset = get_object_or_404(Dataset, id=id)
        stats = {}

        experiments = Experiment.objects.filter(accession_code__in=dataset.data.keys())
        
        # Find all the species for these experiments
        for experiment in experiments:
            species_names = experiment.organisms.values_list('name')
            for species_name in species_names:
                species = stats.get(species_name[0], {"num_experiments": 0, "num_samples": 0})
                species['num_experiments'] = species['num_experiments'] + 1
                stats[species_name[0]] = species

        # Count the samples
        all_sample_accessions = [value[0] for value in dataset.data.values()]
        empty_species = []
        for species in stats.keys():
            samples = Sample.objects.filter(accession_code__in=all_sample_accessions, organism__name=species)
            stats[species]['num_samples'] = len(samples)
            if stats[species]['num_samples'] == 0:
                empty_species.append(species)

        # Delete empty associations
        for species in empty_species:
            del stats[species]

        return Response(stats)

class APITokenView(APIView):
    """
    Return this response to this endpoint with `is_activated: true` to activate this API token.

    You must include an activated token's ID to download processed datasets.
    """

    def get(self, request, id=None):
        """ Create a new token, or fetch a token by its ID. """

        if id:
            token = get_object_or_404(APIToken, id=id)
        else:
            token = APIToken()
            token.save()
        serializer = APITokenSerializer(token)
        return Response(serializer.data)

    def post(self, request, id=None):
        """ Given a token's ID, activate it."""

        id = request.data.get('id', None)
        activated_token = get_object_or_404(APIToken, id=id)
        activated_token.is_activated = request.data.get('is_activated', False)
        activated_token.save()

        serializer = APITokenSerializer(activated_token)
        return Response(serializer.data)

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
        experiments = Experiment.public_objects.filter(**filter_dict)

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
            return Experiment.public_objects.get(pk=pk)
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
    Can also accept a `dataset_id` field instead of a list of accession codes.

    Append the pk or accession_code to the end of this URL to see a detail view.

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        order_by = filter_dict.pop('order_by', None)
        ids = filter_dict.pop('ids', None)
        accession_codes = filter_dict.pop('accession_codes', None)

        if ids is not None:
            ids = [ int(x) for x in ids.split(',')]
            filter_dict['pk__in'] = ids

        if accession_codes is not None:
            accession_codes = accession_codes.split(',')
            filter_dict['accession_code__in'] = accession_codes

        dataset_id = filter_dict.pop('dataset_id', None)
        if dataset_id:
            dataset = get_object_or_404(Dataset, id=dataset_id)
            # Python doesn't provide a prettier way of doing this that I know about.
            filter_dict['accession_code__in'] = [item for sublist in dataset.data.values() for item in sublist]

        samples = Sample.public_objects.filter(**filter_dict)
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
            return Sample.public_objects.get(pk=pk)
        except Sample.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sample = self.get_object(pk)
        serializer = DetailedSampleSerializer(sample)
        return Response(serializer.data)

##
# Processor
##

class ProcessorList(APIView):
    """List all processors."""
    def get(self, request, format=None):
        processors = Processor.objects.all()
        serializer = ProcessorSerializer(processors, many=True)
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
        results = ComputationalResult.public_objects.filter(**filter_dict)

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
        samples = Sample.public_objects.all().values("platform_accession_code", "platform_name").distinct()
        serializer = PlatformSerializer(samples, many=True)
        return Response(serializer.data)

class InstitutionList(APIView):
    """
	Unpaginated list of all the available "institution" information
	"""

    def get(self, request, format=None):
        experiments = Experiment.public_objects.all().values("submitter_institution").distinct()
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

###
# Transcriptome Indices
###

class TranscriptomeIndexDetail(APIView):
    """
    Retrieve the S3 URL and index metadata associated with an OrganismIndex.
    """

    def get(self, request, format=None):
        """
        Gets the S3 url associated with the organism and length, along with other metadata about
        the transcriptome index we have stored. Organism must be specified in underscore-delimited
        uppercase, i.e. "GALLUS_GALLUS". Length must either be "long" or "short"
        """
        params = request.query_params

        # Verify that the required params are present
        errors = dict()
        if "organism" not in params:
            errors["organism"] = "You must specify the organism of the index you want"
        if "length" not in params:
            errors["length"] = "You must specify the length of the transcriptome index"

        if len(errors) > 0:
            raise ValidationError(errors)

        # Get the correct organism index object, serialize it, and return it
        transcription_length = "TRANSCRIPTOME_" + params["length"].upper()
        try:
            organism_index = (OrganismIndex.public_objects.exclude(s3_url__exact="")
                              .distinct("organism__name", "index_type")
                              .get(organism__name=params["organism"],
                                   index_type=transcription_length))
            serializer = OrganismIndexSerializer(organism_index)
            return Response(serializer.data)
        except OrganismIndex.DoesNotExist:
            raise Http404

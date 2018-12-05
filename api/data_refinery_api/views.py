from django.conf import settings
from django.db.models import Count, Prefetch
from django.db.models.aggregates import Avg, Sum
from django.db.models.expressions import F, Q
from django.http import Http404, HttpResponse, HttpResponseRedirect, HttpResponseBadRequest
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
    OrganismIndex,
    ExperimentSampleAssociation,
    OriginalFile,
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
    DatasetDetailsSerializer,
    APITokenSerializer
)

from datetime import timedelta, datetime
from django.utils import timezone

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

    samples__platform_name = \
        django_filters.ModelMultipleChoiceFilter(field_name="samples__platform_name",
                                                 to_field_name="platform_name",
                                                 queryset=Sample.objects.all())
    samples__platform_name.always_filter = False

    class Meta:
        model = Experiment
        fields =    [   'has_publication',
                        'submitter_institution',
                        'technology',
                        'source_first_published',
                        'organisms__name',
                        'samples__platform_name']

    def to_html(self, request, queryset, view):
        # Don't render the FKs in browsable view
        return ''

# Via: https://github.com/encode/django-rest-framework/issues/3905#issuecomment-294391278
class NoMarkupDjangoFilterBackend(DjangoFilterBackend):
    def to_html(self, request, queryset, view):
        # We want this, but currently it incurs a huge performance penality on ChoiceFields with 1000+ choices
        return ''

# ListAPIView is read-only!
class SearchAndFilter(generics.ListAPIView):
    """
    Search and filter for experiments and samples.

    Ex: search/?search=human&has_publication=True

    """

    # Only Experiments with processed objects are exposed
    queryset = Experiment.processed_public_objects.annotate(samples_count=Count('samples')).all()

    # For developing, you can uncomment this to expose everything.
    #queryset = Experiment.objects.annotate(samples_count=Count('samples')).all()

    serializer_class = ExperimentSerializer
    pagination_class = LimitOffsetPagination

    filter_backends = (NoMarkupDjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter)
    filter_class = ExperimentFilter

    # Ordering
    ordering_fields = ('total_samples_count', 'id', 'created_at', 'source_first_published', 'accession_code',)
    samples_count = django_filters.NumberFilter(method='filter_samples_count')
    ordering = ('-total_samples_count',)

    def filter_samples_count(self, queryset, name, value):
        return queryset.filter(total_samples_count=value)

    # via http://www.django-rest-framework.org/api-guide/filtering/#searchfilter
    # '^' Starts-with search.
    # '=' Exact matches.
    # '@' Full-text search.
    # '$' Regex search.
    search_fields = (   'title',
                        '@description',
                        '@accession_code',
                        '@alternate_accession_code',
                        '@protocol_description',
                        '@publication_title',
                        'publication_doi',
                        'publication_authors',
                        'pubmed_id',
                        '@submitter_institution',
                        'experimentannotation__data'
                    )
    filter_fields = ('has_publication', 'platform_name')

    def get_queryset(self):

        # For Prod:
        queryset = Experiment.processed_public_objects.all()

        # For Dev:
        # queryset = Experiment.objects.all()

        # Set up eager loading to avoid N+1 selects
        queryset = self.get_serializer_class().setup_eager_loading(queryset)
        return queryset

    def list(self, request, *args, **kwargs):
        """ Adds counts on certain filter fields to result JSON."""
        response = super(SearchAndFilter, self).list(request, args, kwargs)

        response.data['filters'] = {}
        response.data['filters']['technology'] = {}
        response.data['filters']['publication'] = {}
        response.data['filters']['organism'] = {}
        response.data['filters']['platforms'] = {}

        # Technology
        qs = self.search_queryset(self.get_queryset())
        techs = qs.values('technology').annotate(Count('technology', unique=True))
        for tech in techs:
            if not tech['technology'] or not tech['technology'].strip():
                continue
            response.data['filters']['technology'][tech['technology']] = tech['technology__count']

        # Publication
        pubs = qs.values('has_publication').annotate(Count('has_publication', unique=True))
        for pub in pubs:
            if pub['has_publication']:
                response.data['filters']['publication']['has_publication'] = pub['has_publication__count']

        # Organisms
        organisms = qs.values('organisms__name').annotate(Count('organisms__name', unique=True))
        for organism in organisms:

            # This experiment has no ExperimentOrganism-association, which is bad.
            # This information may still live on the samples though.
            if not organism['organisms__name']:
                continue

            response.data['filters']['organism'][organism['organisms__name']] = organism['organisms__name__count']

        # Platforms
        platforms = qs.values('samples__platform_name').annotate(Count('samples__platform_name', unique=True))
        for plat in platforms:
            if plat['samples__platform_name']:
                response.data['filters']['platforms'][plat['samples__platform_name']] = plat['samples__platform_name__count']

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

    You must also supply `email_address` with `start`, though this will never be serialized back to you.

    """

    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer
    lookup_field = 'id'

    def get_serializer_class(self):
        if 'details' in self.request.query_params:
            return DatasetDetailsSerializer
        return self.serializer_class

    def perform_update(self, serializer):
        """ If `start` is set, fire off the job. Disables dataset data updates after that. """
        old_object = self.get_object()
        old_data = old_object.data
        old_aggregate = old_object.aggregate_by
        already_processing = old_object.is_processing
        new_data = serializer.validated_data

        if old_object.is_processed:
            raise APIException("You may not update Datasets which have already been processed")

        if new_data.get('start'):

            # Make sure we have a valid activated token.
            token_id = self.request.data.get('token_id')
            try:
                token = APIToken.objects.get(id=token_id, is_activated=True)
            except Exception: # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
                raise APIException("You must provide an active API token ID")

            # We could be more aggressive with requirements checking here, but
            # there could be use cases where you don't want to supply an email.
            supplied_email_address = self.request.data.get('email_address', None)

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

                job_sent = False

                obj = serializer.save()
                if supplied_email_address is not None:
                    if obj.email_address != supplied_email_address:
                        obj.email_address = supplied_email_address
                        obj.save()
                try:
                    # Hidden method of non-dispatching for testing purposes.
                    if not self.request.data.get('no_send_job', False):
                        job_sent = send_job(ProcessorPipeline.SMASHER, processor_job)
                    else:
                        # We didn't actually send it, but we also didn't want to.
                        job_sent = True
                except Exception:
                    # job_sent is already false and the exception has
                    # already been logged by send_job, so nothing to
                    # do other than catch the exception.
                    pass

                if not job_sent:
                    raise APIException("Unable to queue download job. Something has gone"
                                       " wrong and we have been notified about it.")

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
        except Exception:
            try:
                return Experiment.public_objects.get(accession_code=pk)
            except Experiment.DoesNotExist:
                raise Http404
            return HttpResponseBadRequest("Bad PK or Accession") 

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

        filter_by = filter_dict.pop('filter_by', None)

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

        if filter_by:
            samples = samples.filter(   Q(sex__contains=filter_by) |
                                        Q(age__contains=filter_by) |
                                        Q(specimen_part__contains=filter_by) |
                                        Q(genotype__contains=filter_by) |
                                        Q(disease__contains=filter_by) |
                                        Q(disease_stage__contains=filter_by) |
                                        Q(cell_line__contains=filter_by) |
                                        Q(treatment__contains=filter_by) |
                                        Q(race__contains=filter_by) |
                                        Q(subject__contains=filter_by) |
                                        Q(compound__contains=filter_by) |
                                        Q(time__contains=filter_by) |
                                        Q(sampleannotation__data__contains=filter_by)
                                    )

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

      Works with required 'limit' and 'offset' params.

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        limit = max(int(filter_dict.pop('limit', 100)), 100)
        offset = int(filter_dict.pop('offset', 0))
        jobs = SurveyJob.objects.filter(**filter_dict)[offset:(offset + limit)]
        serializer = SurveyJobSerializer(jobs, many=True)
        return Response(serializer.data)

class DownloaderJobList(PaginatedAPIView):
    """
    List of all DownloaderJob
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        limit = max(int(filter_dict.pop('limit', 100)), 100)
        offset = int(filter_dict.pop('offset', 0))
        jobs = DownloaderJob.objects.filter(**filter_dict)[offset: offset + limit]
        serializer = DownloaderJobSerializer(jobs, many=True)
        return Response(serializer.data)

class ProcessorJobList(PaginatedAPIView):
    """
    List of all ProcessorJobs
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        limit = max(int(filter_dict.pop('limit', 100)), 100)
        offset = int(filter_dict.pop('offset', 0))
        jobs = ProcessorJob.objects.filter(**filter_dict)[offset: offset + limit]
        serializer = ProcessorJobSerializer(jobs, many=True)
        return Response(serializer.data)

###
# Statistics
###

class Stats(APIView):
    """
    Statistics about the health of the system.

    ?range=week  includes statics for the last week
    """

    def get(self, request, format=None):
        range_param = request.query_params.dict().pop('range', None)

        data = {}
        data['survey_jobs'] = self._get_job_stats(SurveyJob.objects, range_param)
        data['downloader_jobs'] = self._get_job_stats(DownloaderJob.objects, range_param)
        data['processor_jobs'] = self._get_job_stats(ProcessorJob.objects, range_param)
        data['samples'] = self._get_object_stats(Sample.objects, range_param)
        data['experiments'] = self._get_object_stats(Experiment.objects, range_param)
        data['total_downloaded'] = self._get_total_downloaded()

        return Response(data)

    def _get_total_downloaded(self):

        total_size = OriginalFile.objects.filter(is_downloaded=True).aggregate(Sum('size_in_bytes'))
        return total_size['size_in_bytes__sum']

    def _get_job_stats(self, jobs, range_param):
        result = {
            'total': jobs.count(),
            'pending': jobs.filter(start_time__isnull=True).count(),
            'completed': jobs.filter(end_time__isnull=False).count(),
            'open': jobs.filter(start_time__isnull=False, end_time__isnull=True, success__isnull=True).count(),
            # via https://stackoverflow.com/questions/32520655/get-average-of-difference-of-datetime-fields-in-django
            'average_time': jobs.filter(start_time__isnull=False, end_time__isnull=False).aggregate(
                average_time=Avg(F('end_time') - F('start_time')))['average_time']
        }

        if not result['average_time']:
            result['average_time'] = 0
        else:
            result['average_time'] = result['average_time'].total_seconds()

        if range_param:
            result['timeline'] = self._jobs_timeline(jobs, range_param)

        return result

    def _get_object_stats(self, objects, range_param):
        result = {
            'total': objects.count()
        }

        if range_param:
            result['timeline'] = self._created_timeline(objects, range_param)

        return result

    def _get_time_intervals(self, range_param):
        interval_timedelta = {
            'day': timedelta(days=1),
            'week': timedelta(weeks=1),
            'month': timedelta(weeks=4),
            'year': timedelta(weeks=52)
        }
        interval_timestep = {
            'day': timedelta(hours=1),
            'week': timedelta(days=1),
            'month': timedelta(days=2),
            'year': timedelta(weeks=4)
        }

        current_date = datetime.now(tz=timezone.utc)
        time_step = interval_timestep.get(range_param)
        start_date = current_date - interval_timedelta.get(range_param)

        intervals = [(current_date - time_step*(i+1), current_date - time_step*i)
                     for i in range(100) if current_date - time_step*(i+1) > start_date]
        return intervals[::-1]

    def _get_job_interval(self, jobs, start, end):
        filtered_jobs = jobs.filter(created_at__gte=start, created_at__lte=end)
        pending = filtered_jobs and jobs.filter(start_time__isnull=True)
        failed = filtered_jobs and jobs.filter(success=False)
        completed = filtered_jobs and jobs.filter(success=True)
        open = filtered_jobs and jobs.filter(success__isnull=True)

        return {
            'start': start,
            'end': end,
            'total': filtered_jobs.count(),
            'completed': completed.count(),
            'pending': pending.count(),
            'failed': failed.count(),
            'open': open.count()
        }

    def _jobs_timeline(self, jobs, range_param):
        return [self._get_job_interval(jobs, start, end) for (start, end) in self._get_time_intervals(range_param)]

    def _created_timeline(self, objects, range_param):
        results = []
        for start, end in self._get_time_intervals(range_param):
            total = objects.filter(created_at__gte=start, created_at__lte=end).count()
            stats = {
                'start': start, 
                'end': end, 
                'total': total
            }
            results.append(stats)
        return results

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
                              .get(organism__name=params["organism"].upper(),
                                   index_type=transcription_length))
            serializer = OrganismIndexSerializer(organism_index)
            return Response(serializer.data)
        except OrganismIndex.DoesNotExist:
            raise Http404

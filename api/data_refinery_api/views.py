from datetime import timedelta, datetime
import requests
import nomad
from typing import Dict
from itertools import groupby
from re import match
from django.conf import settings
from django.views.decorators.cache import cache_page
from django.db.models import Count, Prefetch, DateTimeField
from django.db.models.functions import Trunc
from django.db.models.aggregates import Avg, Sum
from django.db.models.expressions import F, Q
from django.http import Http404, HttpResponse, HttpResponseRedirect, HttpResponseBadRequest
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django_elasticsearch_dsl_drf.constants import (
    LOOKUP_FILTER_TERMS,
    LOOKUP_FILTER_RANGE,
    LOOKUP_FILTER_PREFIX,
    LOOKUP_FILTER_WILDCARD,
    LOOKUP_QUERY_IN,
    LOOKUP_QUERY_GT,
    LOOKUP_QUERY_GTE,
    LOOKUP_QUERY_LT,
    LOOKUP_QUERY_LTE,
    LOOKUP_QUERY_EXCLUDE,
)
from django_elasticsearch_dsl_drf.viewsets import DocumentViewSet
from django_elasticsearch_dsl_drf.filter_backends import (
    FilteringFilterBackend,
    IdsFilterBackend,
    OrderingFilterBackend,
    DefaultOrderingFilterBackend,
    CompoundSearchFilterBackend,
    FacetedSearchFilterBackend
)
from django_filters.rest_framework import DjangoFilterBackend
import django_filters
from elasticsearch_dsl import TermsFacet, DateHistogramFacet
from rest_framework import status, filters, generics, mixins
from rest_framework.exceptions import APIException, NotFound
from rest_framework.exceptions import ValidationError
from rest_framework.pagination import LimitOffsetPagination
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.settings import api_settings
from rest_framework.views import APIView

from data_refinery_api.serializers import (
    ComputationalResultSerializer,
    ComputationalResultWithUrlSerializer,
    DetailedExperimentSerializer,
    DetailedSampleSerializer,
    ExperimentSerializer,
    InstitutionSerializer,
    OrganismIndexSerializer,
    OrganismSerializer,
    PlatformSerializer,
    ProcessorSerializer,
    SampleSerializer,
    CompendiaSerializer,
    CompendiaWithUrlSerializer,
    QNTargetSerializer,
    ComputedFileListSerializer,

    # Job
    DownloaderJobSerializer,
    ProcessorJobSerializer,
    SurveyJobSerializer,

    # Dataset
    APITokenSerializer,
    CreateDatasetSerializer,
    DatasetSerializer,
)
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    APIToken,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    DownloaderJob,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    OriginalFile,
    Processor,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SurveyJob,
)
from data_refinery_common.models.documents import (
    ExperimentDocument
)
from data_refinery_common.utils import (
    get_env_variable,
    get_active_volumes,
    get_nomad_jobs_breakdown,
)
from data_refinery_common.logging import get_and_configure_logger
from .serializers import ExperimentDocumentSerializer

from django.utils.decorators import method_decorator
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

logger = get_and_configure_logger(__name__)

##
# Variables
##

JOB_CREATED_AT_CUTOFF = datetime(2019, 6, 5, tzinfo=timezone.utc)

##
# ElasticSearch
##
from django_elasticsearch_dsl_drf.pagination import LimitOffsetPagination as ESLimitOffsetPagination
from six import iteritems

class FacetedSearchFilterBackendExtended(FacetedSearchFilterBackend):
    def aggregate(self, request, queryset, view):
        """Extends FacetedSearchFilterBackend to add additional metrics to each bucket
        https://github.com/barseghyanartur/django-elasticsearch-dsl-drf/blob/master/src/django_elasticsearch_dsl_drf/filter_backends/faceted_search.py#L19

        We have the downloadable sample accession codes indexed for each experiment.
        The cardinality metric, returns the number of unique samples for each bucket.
        However it's just an approximate
        https://www.elastic.co/guide/en/elasticsearch/reference/current/search-aggregations-metrics-cardinality-aggregation.html#_counts_are_approximate
        I used the highest possible precision threshold, but this might increase the amount
        of memory used.
        """
        facets = self.construct_facets(request, view)
        for field, facet in iteritems(facets):
            agg = facet['facet'].get_aggregation()
            queryset.aggs.bucket(field, agg)\
                .metric('total_samples', 'cardinality', field='downloadable_samples', precision_threshold=40000)
        return queryset


##
# ElasticSearch powered Search and Filter
##
@method_decorator(name='list', decorator=swagger_auto_schema(
manual_parameters=[
    openapi.Parameter(
        name='technology', in_=openapi.IN_QUERY, type=openapi.TYPE_STRING,
        description="Allows filtering the results by technology, can have multiple values. Eg: `?technology=microarray&technology=rna-seq`",
    ),
    openapi.Parameter(
        name='has_publication', in_=openapi.IN_QUERY, type=openapi.TYPE_STRING,
        description="Filter the results that have associated publications with `?has_publication=true`",
    ),
    openapi.Parameter(
        name='platform', in_=openapi.IN_QUERY,
        type=openapi.TYPE_STRING,
        description="Allows filtering the results by platform, this parameter can have multiple values.",
    ),
    openapi.Parameter(
        name='organism', in_=openapi.IN_QUERY,
        type=openapi.TYPE_STRING,
        description="Allows filtering the results by organism, this parameter can have multiple values.",
    ),
    openapi.Parameter(
        name='num_processed_samples', in_=openapi.IN_QUERY,
        type=openapi.TYPE_NUMBER,
        description="Use ElasticSearch queries to specify the number of processed samples of the results",
    ),
],
operation_description="""
Use this endpoint to search among the experiments.

This is powered by ElasticSearch, information regarding advanced usages of the
filters can be found in the [Django-ES-DSL-DRF docs](https://django-elasticsearch-dsl-drf.readthedocs.io/en/0.17.1/filtering_usage_examples.html#filtering)

There's an additional field in the response named `facets` that contain stats on the number of results per filter type.

Example Requests:
```
?search=medulloblastoma
?id=1
?search=medulloblastoma&technology=microarray&has_publication=true
?ordering=source_first_published
```
"""))
class ExperimentDocumentView(DocumentViewSet):
    """ ElasticSearch powered experiment search. """
    document = ExperimentDocument
    serializer_class = ExperimentDocumentSerializer
    pagination_class = ESLimitOffsetPagination

    # Filter backends provide different functionality we want
    filter_backends = [
        FilteringFilterBackend,
        OrderingFilterBackend,
        DefaultOrderingFilterBackend,
        CompoundSearchFilterBackend,
        FacetedSearchFilterBackendExtended
    ]

    # Primitive
    lookup_field = 'id'

    # Define search fields
    # Is this exhaustive enough?
    search_fields = {
        'title': {'boost': 10},
        'publication_authors': {'boost': 8},  # "People will search themselves"
        'publication_title':  {'boost': 5},
        'submitter_institution': {'boost': 3},
        'description':  {'boost': 2},
        'accession_code': None,
        'alternate_accession_code': None,
        'publication_doi': None,
        'pubmed_id': None,
        'sample_metadata_fields': None,
        'platform_names': None
    }

    # Define filtering fields
    filter_fields = {
        'id': {
            'field': '_id',
            'lookups': [
                LOOKUP_FILTER_RANGE,
                LOOKUP_QUERY_IN
            ],
        },
        'technology': 'technology',
        'has_publication': 'has_publication',
        'platform': 'platform_accession_codes',
        'organism': 'organism_names',
        'num_processed_samples': {
            'field': 'num_processed_samples',
            'lookups': [
                LOOKUP_FILTER_RANGE,
                LOOKUP_QUERY_IN,
                LOOKUP_QUERY_GT
            ],
        },
        'num_downloadable_samples': {
            'field': 'num_downloadable_samples',
            'lookups': [
                LOOKUP_FILTER_RANGE,
                LOOKUP_QUERY_IN,
                LOOKUP_QUERY_GT
            ],
        }
    }

    # Define ordering fields
    ordering_fields = {
        'id': 'id',
        'title': 'title.raw',
        'description': 'description.raw',
        'num_total_samples': 'num_total_samples',
        'num_downloadable_samples': 'num_downloadable_samples',
        'source_first_published': 'source_first_published'
    }

    # Specify default ordering
    ordering = ('_score', '-num_total_samples', 'id', 'title', 'description', '-source_first_published')

    # Facets (aka Aggregations) provide statistics about the query result set in the API response.
    # More information here: https://github.com/barseghyanartur/django-elasticsearch-dsl-drf/blob/03a3aa716db31868ca3a71340513a993741a4177/src/django_elasticsearch_dsl_drf/filter_backends/faceted_search.py#L24
    faceted_search_fields = {
        'technology': {
            'field': 'technology',
            'facet': TermsFacet,
            'enabled': True  # These are enabled by default, which is more expensive but more simple.
        },
        'organism_names': {
            'field': 'organism_names',
            'facet': TermsFacet,
            'enabled': True,
            'options': {
                'size': 999999
            }
        },
        'platform_accession_codes': {
            'field': 'platform_accession_codes',
            'facet': TermsFacet,
            'enabled': True,
            'global': False,
            'options': {
                'size': 999999
            }
        },
        'has_publication': {
            'field': 'has_publication',
            'facet': TermsFacet,
            'enabled': True,
            'global': False,
        },
        # We don't actually need any "globals" to drive our web frontend,
        # but we'll leave them available but not enabled by default, as they're
        # expensive.
        'technology_global': {
            'field': 'technology',
            'facet': TermsFacet,
            'enabled': False,
            'global': True
        },
        'organism_names_global': {
            'field': 'organism_names',
            'facet': TermsFacet,
            'enabled': False,
            'global': True,
            'options': {
                'size': 999999
            }
        },
        'platform_names_global': {
            'field': 'platform_names',
            'facet': TermsFacet,
            'enabled': False,
            'global': True,
            'options': {
                'size': 999999
            }
        },
        'has_publication_global': {
            'field': 'platform_names',
            'facet': TermsFacet,
            'enabled': False,
            'global': True,
        },
    }
    faceted_search_param = 'facet'

    def list(self, request, *args, **kwargs):
        response = super(ExperimentDocumentView, self).list(request, args, kwargs)
        response.data['facets'] = self.transform_es_facets(response.data['facets'])
        return response

    def transform_es_facets(self, facets):
        """Transforms Elastic Search facets into a set of objects where each one corresponds
        to a filter group. Example:

        { technology: {rna-seq: 254, microarray: 8846, unknown: 0} }

        Which means the users could attach `?technology=rna-seq` to the url and expect 254
        samples returned in the results.
        """
        result = {}
        for field, facet in iteritems(facets):
            filter_group = {}
            for bucket in facet['buckets']:
                if field == 'has_publication':
                    filter_group[bucket['key_as_string']] = bucket['total_samples']['value']
                else:
                    filter_group[bucket['key']] = bucket['total_samples']['value']
            result[field] = filter_group
        return result

##
# Dataset
##

class CreateDatasetView(generics.CreateAPIView):
    """ Creates and returns new Datasets. """
    queryset = Dataset.objects.all()
    serializer_class = CreateDatasetSerializer

@method_decorator(name='get', decorator=swagger_auto_schema(operation_description="View a single Dataset.",manual_parameters=[
openapi.Parameter(
    name='details', in_=openapi.IN_QUERY, type=openapi.TYPE_BOOLEAN,
    description="When set to `True`, additional fields will be included in the response with details about the experiments in the dataset. This is used mostly on the dataset page in www.refine.bio",
)]))
@method_decorator(name='patch', decorator=swagger_auto_schema(auto_schema=None)) # partial updates not supported
@method_decorator(name='put', decorator=swagger_auto_schema(operation_description="""
Modify an existing Dataset.

Set `start` to `true` along with a valid activated API token (from `/token/`) to begin smashing and delivery.

You must also supply `email_address` with `start`, though this will never be serialized back to you.
"""))
class DatasetView(generics.RetrieveUpdateAPIView):
    """ View and modify a single Dataset. """
    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer
    lookup_field = 'id'

    @staticmethod
    def _should_display_on_engagement_bot(email: str) -> bool:
        return email is not None \
            and email.find("cansav09") != 0 \
            and email.find("arielsvn") != 0 \
            and email.find("jaclyn.n.taroni") != 0 \
            and email.find("kurt.wheeler") != 0 \
            and email.find("greenescientist") != 0 \
            and email.find("@alexslemonade.org") == -1 \
            and email.find("miserlou") != 0 \
            and email.find("michael.zietz@gmail.com") != 0 \
            and email.find("d.prasad") != 0 \
            and email.find("daniel.himmelstein@gmail.com") != 0 \
            and email.find("dv.prasad991@gmail.com") != 0

    def get_serializer_context(self):
        """
        Extra context provided to the serializer class.
        """
        serializer_context = super(DatasetView, self).get_serializer_context()
        token_id = self.request.META.get('HTTP_API_KEY', None)
        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return {**serializer_context, 'token': token}
        except Exception:  # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return serializer_context

    def perform_update(self, serializer):
        """ If `start` is set, fire off the job. Disables dataset data updates after that. """
        old_object = self.get_object()
        old_data = old_object.data
        old_aggregate = old_object.aggregate_by
        already_processing = old_object.is_processing
        new_data = serializer.validated_data

        qn_organisms = Organism.get_objects_with_qn_targets()

        # We convert 'ALL' into the actual accession codes given
        for key in new_data['data'].keys():
            accessions = new_data['data'][key]
            if accessions == ["ALL"]:
                experiment = get_object_or_404(Experiment, accession_code=key)

                sample_codes = list(experiment.samples.filter(is_processed=True, organism__in=qn_organisms).values_list('accession_code', flat=True))
                new_data['data'][key] = sample_codes

        if old_object.is_processed:
            raise APIException("You may not update Datasets which have already been processed")

        if new_data.get('start'):

            # Make sure we have a valid activated token.
            token_id = self.request.data.get('token_id', None)

            if not token_id:
                token_id = self.request.META.get('HTTP_API_KEY', None)

            try:
                token = APIToken.objects.get(id=token_id, is_activated=True)
            except Exception: # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
                raise APIException("You must provide an active API token ID")

            supplied_email_address = self.request.data.get('email_address', None)
            email_ccdl_ok = self.request.data.get('email_ccdl_ok', False)

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
                if email_ccdl_ok:
                    obj.email_ccdl_ok = email_ccdl_ok
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

                if settings.RUNNING_IN_CLOUD and settings.ENGAGEMENTBOT_WEBHOOK is not None \
                   and DatasetView._should_display_on_engagement_bot(supplied_email_address):
                    try:
                        try:
                            remote_ip = get_client_ip(self.request)
                            city = requests.get('https://ipapi.co/' + remote_ip + '/json/', timeout=10).json()['city']
                        except Exception:
                            city = "COULD_NOT_DETERMINE"

                        new_user_text = "New user " + supplied_email_address + " from " + city + " [" + remote_ip + "] downloaded a dataset! (" + str(old_object.id) + ")"
                        webhook_url = settings.ENGAGEMENTBOT_WEBHOOK
                        slack_json = {
                            "channel": "ccdl-general", # Move to robots when we get sick of these
                            "username": "EngagementBot",
                            "icon_emoji": ":halal:",
                            "attachments":[
                                {   "color": "good",
                                    "text": new_user_text
                                }
                            ]
                        }
                        response = requests.post(
                            webhook_url,
                            json=slack_json,
                            headers={'Content-Type': 'application/json'},
                            timeout=10
                        )
                    except Exception as e:
                        # It doens't really matter if this didn't work
                        logger.error(e)
                        pass

                return obj

        # Don't allow critical data updates to jobs that have already been submitted,
        # but do allow email address updating.
        if already_processing:
            serializer.validated_data['data'] = old_data
            serializer.validated_data['aggregate_by'] = old_aggregate
        serializer.save()

class CreateApiTokenView(generics.CreateAPIView):
    """
    token_create

    There're several endpoints like [/dataset](#tag/dataset) and [/results](#tag/results) that return
    S3 urls where users can download the files we produce, however in order to get those files people
    need to accept our terms of use by creating a token and activating it.

    ```
    POST /token
    PUT /token/{token-id} is_active=True
    ```

    The token id needs to be sent on the `API_KEY` header on http requests.

    References
    - [https://github.com/AlexsLemonade/refinebio/issues/731]()
    - [https://github.com/AlexsLemonade/refinebio-frontend/issues/560]()
    """
    model = APIToken
    serializer_class = APITokenSerializer

@method_decorator(name='patch', decorator=swagger_auto_schema(auto_schema=None))
class APITokenView(generics.RetrieveUpdateAPIView):
    """
    Read and modify Api Tokens.

    get:
    Return details about a specific token.

    put:
    This can be used to activate a specific token by sending `is_activated: true`.
    """
    model = APIToken
    lookup_field = 'id'
    queryset = APIToken.objects.all()
    serializer_class = APITokenSerializer

##
# Experiments
##

class ExperimentList(generics.ListAPIView):
    """ Paginated list of all experiments. Advanced filtering can be done with the `/search` endpoint. """
    model = Experiment
    queryset = Experiment.public_objects.all()
    serializer_class = ExperimentSerializer
    filter_backends = (DjangoFilterBackend,)
    filterset_fields = (
        'title',
        'description',
        'accession_code',
        'alternate_accession_code',
        'source_database',
        'source_url',
        'has_publication',
        'publication_title',
        'publication_doi',
        'pubmed_id',
        'organisms',
        'submitter_institution',
        'created_at',
        'last_modified',
        'source_first_published',
        'source_last_modified',
    )

class ExperimentDetail(generics.RetrieveAPIView):
    """ Retrieve details for an experiment given it's accession code """
    lookup_field = "accession_code"
    queryset = Experiment.public_objects.all()
    serializer_class = DetailedExperimentSerializer

##
# Samples
##

@method_decorator(name='get', decorator=swagger_auto_schema(manual_parameters=[
    openapi.Parameter(
        name='dataset_id', in_=openapi.IN_QUERY,
        type=openapi.TYPE_STRING,
        description="Filters the result and only returns samples that are added to a dataset.",
    ),
    openapi.Parameter(
        name='experiment_accession_code', in_=openapi.IN_QUERY,
        type=openapi.TYPE_STRING,
        description="Filters the result and only returns only the samples associated with an experiment accession code.",
    ),
    openapi.Parameter(
        name='accession_codes', in_=openapi.IN_QUERY,
        type=openapi.TYPE_STRING,
        description="Provide a list of sample accession codes separated by commas and the endpoint will only return information about these samples.",
    ),
]))
class SampleList(generics.ListAPIView):
    """ Returns detailed information about Samples """
    model = Sample
    serializer_class = DetailedSampleSerializer
    filter_backends = (filters.OrderingFilter, DjangoFilterBackend)
    ordering_fields = '__all__'
    ordering = ('-is_processed')
    filterset_fields = (
        'title',
        'organism',
        'source_database',
        'source_archive_url',
        'has_raw',
        'platform_name',
        'technology',
        'manufacturer',
        'sex',
        'age',
        'specimen_part',
        'genotype',
        'disease',
        'disease_stage',
        'cell_line',
        'treatment',
        'race',
        'subject',
        'compound',
        'time',
        'is_processed',
        'is_public'
    )

    def get_queryset(self):
        """
        ref https://www.django-rest-framework.org/api-guide/filtering/#filtering-against-query-parameters
        """
        queryset = Sample.public_objects \
            .prefetch_related('organism') \
            .prefetch_related('results') \
            .prefetch_related('results__processor') \
            .prefetch_related('results__computationalresultannotation_set') \
            .prefetch_related('results__computedfile_set') \
            .filter(**self.get_query_params_filters())

        # case insensitive search https://docs.djangoproject.com/en/2.1/ref/models/querysets/#icontains
        filter_by = self.request.query_params.get('filter_by', None)
        if filter_by:
            queryset = queryset.filter(Q(title__icontains=filter_by) |
                                       Q(sex__icontains=filter_by) |
                                       Q(age__icontains=filter_by) |
                                       Q(specimen_part__icontains=filter_by) |
                                       Q(genotype__icontains=filter_by) |
                                       Q(disease__icontains=filter_by) |
                                       Q(disease_stage__icontains=filter_by) |
                                       Q(cell_line__icontains=filter_by) |
                                       Q(treatment__icontains=filter_by) |
                                       Q(race__icontains=filter_by) |
                                       Q(subject__icontains=filter_by) |
                                       Q(compound__icontains=filter_by) |
                                       Q(time__icontains=filter_by))

        return queryset

    def get_query_params_filters(self):
        """ We do advanced filtering on the queryset depending on the query parameters.
            This returns the parameters that should be used for that. """
        filter_dict = dict()

        ids = self.request.query_params.get('ids', None)
        if ids is not None:
            ids = [int(x) for x in ids.split(',')]
            filter_dict['pk__in'] = ids

        experiment_accession_code = self.request.query_params.get('experiment_accession_code', None)
        if experiment_accession_code:
            experiment = get_object_or_404(Experiment.objects.values('id'), accession_code=experiment_accession_code)
            filter_dict['experiments__in'] = [experiment['id']]

        accession_codes = self.request.query_params.get('accession_codes', None)
        if accession_codes:
            accession_codes = accession_codes.split(',')
            filter_dict['accession_code__in'] = accession_codes

        dataset_id = self.request.query_params.get('dataset_id', None)
        if dataset_id:
            dataset = get_object_or_404(Dataset, id=dataset_id)
            # Python doesn't provide a prettier way of doing this that I know about.
            filter_dict['accession_code__in'] = [item for sublist in dataset.data.values() for item in sublist]

        # Accept Organism in both name and ID form
        organism_name = self.request.query_params.get('organism__name', None)
        if organism_name:
            filter_dict['organism__name'] = organism_name

        return filter_dict

class SampleDetail(generics.RetrieveAPIView):
    """ Retrieve the details for a Sample given it's accession code """
    lookup_field = "accession_code"
    queryset = Sample.public_objects.all()
    serializer_class = DetailedSampleSerializer

##
# Processor
##

class ProcessorList(generics.ListAPIView):
    """List all processors."""
    queryset = Processor.objects.all()
    serializer_class = ProcessorSerializer


##
# Results
##

class ComputationalResultsList(generics.ListAPIView):
    """
    computational_results_list

    This lists all `ComputationalResult`. Each one contains meta-information about the output of a computer process. (Ex Salmon).

    This can return valid S3 urls if a valid [token](#tag/token) is sent in the header `HTTP_API_KEY`.
    """
    queryset = ComputationalResult.public_objects.all()

    def get_serializer_class(self):
        token_id = self.request.META.get('HTTP_API_KEY', None)

        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return ComputationalResultWithUrlSerializer
        except Exception: # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return ComputationalResultSerializer

    def filter_queryset(self, queryset):
        filter_dict = self.request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        return queryset.filter(**filter_dict)

##
# Search Filter Models
##

class OrganismList(generics.ListAPIView):
    """
	Unpaginated list of all the available organisms.
	"""
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer
    paginator = None

class PlatformList(generics.ListAPIView):
    """
	Unpaginated list of all the available "platform" information
	"""
    serializer_class = PlatformSerializer
    paginator = None

    def get_queryset(self):
        return Sample.public_objects.all().values("platform_accession_code", "platform_name").distinct()

class InstitutionList(generics.ListAPIView):
    """
	Unpaginated list of all the available "institution" information
	"""
    serializer_class = InstitutionSerializer
    paginator = None

    def get_queryset(self):
        return Experiment.public_objects.all().values("submitter_institution").distinct()

##
# Jobs
##

class SurveyJobList(generics.ListAPIView):
    """
    List of all SurveyJob.
    """
    model = SurveyJob
    queryset = SurveyJob.objects.all()
    serializer_class = SurveyJobSerializer
    filter_backends = (DjangoFilterBackend, filters.OrderingFilter,)
    filterset_fields = SurveyJobSerializer.Meta.fields
    ordering_fields = ('id', 'created_at')
    ordering = ('-id',)

class DownloaderJobList(generics.ListAPIView):
    """
    List of all DownloaderJob
    """
    model = DownloaderJob
    queryset = DownloaderJob.objects.all()
    serializer_class = DownloaderJobSerializer
    filter_backends = (DjangoFilterBackend, filters.OrderingFilter,)
    filterset_fields = DownloaderJobSerializer.Meta.fields
    ordering_fields = ('id', 'created_at')
    ordering = ('-id',)

class ProcessorJobList(generics.ListAPIView):
    """
    List of all ProcessorJobs.
    """
    model = ProcessorJob
    queryset = ProcessorJob.objects.all()
    serializer_class = ProcessorJobSerializer
    filter_backends = (DjangoFilterBackend, filters.OrderingFilter,)
    filterset_fields = ProcessorJobSerializer.Meta.fields
    ordering_fields = ('id', 'created_at')
    ordering = ('-id',)

###
# Statistics
###

class AboutStats(APIView):
    """ Returns general stats for the site, used in the about page """
    @method_decorator(cache_page(10 * 60))
    def get(self, request, version, format=None):
        # static values for now
        dummy = request.query_params.dict().pop('dummy', None)
        if dummy:
            # add a dummy response, calculated these on 09/25/2019
            result = {
                'samples_available': 904953 + 391022,
                'total_size_in_bytes': 832195361132962,
                'supported_organisms': 43 + 159,
                'experiments_processed': 35785 + 8661
            }
            return Response(result)

        result = {
            'samples_available': self._get_samples_available(),
            'total_size_in_bytes': OriginalFile.objects.aggregate(total_size=Sum('size_in_bytes'))['total_size'],
            'supported_organisms': self._get_supported_organisms(),
            'experiments_processed': self._get_experiments_processed()
        }
        return Response(result)

    def _get_experiments_processed(self):
        """ total experiments with at least one sample processed """
        experiments_with_sample_processed = Experiment.objects\
            .annotate(
                processed_samples_count=Count('samples', filter=Q(samples__is_processed=True)),
            )\
            .filter(Q(processed_samples_count__gt=1))\
            .count()
        experiments_with_sample_quant = ComputedFile.objects\
            .filter(filename='quant.sf', result__samples__is_processed=False)\
            .values_list('result__samples__experiments', flat=True)\
            .distinct()\
            .count()
        return experiments_with_sample_processed + experiments_with_sample_quant

    def _get_supported_organisms(self):
        """ count organisms with qn targets or that have at least one sample with quant files """
        organisms_with_qn_targets = Organism.objects.filter(qn_target__isnull=False).count()
        organisms_without_qn_targets = Organism.objects.filter(
            qn_target__isnull=True,
            sample__is_processed=False,
            sample__technology='RNA-SEQ',
            sample__results__computedfile__filename='quant.sf'
        ).distinct().count()
        return organisms_with_qn_targets + organisms_without_qn_targets

    def _get_samples_available(self):
        """ count the total number of samples that are processed or that have a quant.sf file associated with them """
        processed_samples = Sample.objects.filter(is_processed=True).count()
        unprocessed_samples_with_quant = Sample.objects.filter(
            is_processed=False,
            technology='RNA-SEQ',
            results__computedfile__filename='quant.sf'
        ).distinct().count()
        return processed_samples + unprocessed_samples_with_quant

class Stats(APIView):
    """ Statistics about the health of the system. """

    @swagger_auto_schema(manual_parameters=[openapi.Parameter(
        name='range', in_=openapi.IN_QUERY, type=openapi.TYPE_STRING,
        description="Specify a range from which to calculate the possible options",
        enum=('day', 'week', 'month', 'year',)
    )])
    @method_decorator(cache_page(10 * 60))
    def get(self, request, version, format=None):
        range_param = request.query_params.dict().pop('range', None)
        cached_stats = Stats.calculate_stats(range_param)
        return Response(cached_stats)

    @classmethod
    def calculate_stats(cls, range_param):
        data = {}
        data['generated_on'] = timezone.now()
        data['survey_jobs'] = cls._get_job_stats(SurveyJob.objects, range_param)
        data['downloader_jobs'] = cls._get_job_stats(DownloaderJob.objects, range_param)
        data['processor_jobs'] = cls._get_job_stats(ProcessorJob.objects, range_param)
        data['experiments'] = cls._get_object_stats(Experiment.objects, range_param)

        # processed and unprocessed samples stats
        data['unprocessed_samples'] = cls._get_object_stats(Sample.objects.filter(is_processed=False), range_param, 'last_modified')
        data['processed_samples'] = cls._get_object_stats(Sample.processed_objects, range_param, 'last_modified')
        data['processed_samples']['last_hour'] = cls._samples_processed_last_hour()

        data['processed_samples']['technology'] = {}
        techs = Sample.processed_objects.values('technology').annotate(count=Count('technology'))
        for tech in techs:
            if not tech['technology'] or not tech['technology'].strip():
                continue
            data['processed_samples']['technology'][tech['technology']] = tech['count']

        data['processed_samples']['organism'] = {}
        organisms = Sample.processed_objects.values('organism__name').annotate(count=Count('organism__name'))
        for organism in organisms:
            if not organism['organism__name']:
                continue
            data['processed_samples']['organism'][organism['organism__name']] = organism['count']

        data['processed_experiments'] = cls._get_object_stats(Experiment.processed_public_objects)
        data['active_volumes'] = list(get_active_volumes())
        data['dataset'] = cls._get_dataset_stats(range_param)

        if range_param:
            data['input_data_size'] = cls._get_input_data_size()
            data['output_data_size'] = cls._get_output_data_size()

        data.update(get_nomad_jobs_breakdown())

        return data

    EMAIL_USERNAME_BLACKLIST = [
        'arielsvn',
        'cansav09',
        'd.prasad',
        'daniel.himmelstein',
        'dv.prasad991',
        'greenescientist',
        'jaclyn.n.taroni',
        'kurt.wheeler91',
        'michael.zietz',
        'miserlou'
    ]

    @classmethod
    def _get_dataset_stats(cls, range_param):
        """Returns stats for processed datasets"""
        filter_query = Q()
        for username in Stats.EMAIL_USERNAME_BLACKLIST:
            filter_query = filter_query | Q(email_address__startswith=username)
        filter_query = filter_query | Q(email_address__endswith='@alexslemonade.org')
        processed_datasets = Dataset.objects.filter(is_processed=True, email_address__isnull=False).exclude(filter_query)
        result = processed_datasets.aggregate(
            total=Count('id'),
            aggregated_by_experiment=Count('id', filter=Q(aggregate_by='EXPERIMENT')),
            aggregated_by_species=Count('id', filter=Q(aggregate_by='SPECIES')),
            scale_by_none=Count('id', filter=Q(scale_by='NONE')),
            scale_by_minmax=Count('id', filter=Q(scale_by='MINMAX')),
            scale_by_standard=Count('id', filter=Q(scale_by='STANDARD')),
            scale_by_robust=Count('id', filter=Q(scale_by='ROBUST')),
        )

        if range_param:
            # We don't save the dates when datasets are processed, but we can use
            # `last_modified`, since datasets aren't modified again after they are processed
            result['timeline'] = cls._get_intervals(
                processed_datasets,
                range_param,
                'last_modified'
            ).annotate(
                total=Count('id'),
                total_size=Sum('size_in_bytes')
            )
        return result

    @classmethod
    def _samples_processed_last_hour(cls):
        current_date = datetime.now(tz=timezone.utc)
        start = current_date - timedelta(hours=1)
        return Sample.processed_objects.filter(last_modified__range=(start, current_date)).count()

    @classmethod
    def _get_input_data_size(cls):
        total_size = OriginalFile.objects.filter(
            sample__is_processed=True  # <-- SLOW
        ).aggregate(
            Sum('size_in_bytes')
        )
        return total_size['size_in_bytes__sum'] if total_size['size_in_bytes__sum'] else 0

    @classmethod
    def _get_output_data_size(cls):
        total_size = ComputedFile.public_objects.all().filter(
            s3_bucket__isnull=False,
            s3_key__isnull=True
        ).aggregate(
            Sum('size_in_bytes')
        )
        return total_size['size_in_bytes__sum'] if total_size['size_in_bytes__sum'] else 0

    @classmethod
    def _get_job_stats(cls, jobs, range_param):
        start_filter = Q()

        if range_param:
            current_date = datetime.now(tz=timezone.utc)
            start_date = {
                'day': current_date - timedelta(days=1),
                'week': current_date - timedelta(weeks=1),
                'month': current_date - timedelta(days=30),
                'year': current_date - timedelta(days=365)
            }.get(range_param)
            start_filter = start_filter | Q(start_time__gte=start_date) | Q(start_time__isnull=True)

        result = jobs.filter(start_filter).aggregate(
            total=Count('id'),
            successful=Count('id', filter=Q(success=True)),
            failed=Count('id', filter=Q(success=False)),
            pending=Count('id', filter=Q(start_time__isnull=True,
                                         success__isnull=True,
                                         created_at__gt=JOB_CREATED_AT_CUTOFF)),
            open=Count('id', filter=Q(start_time__isnull=False,
                                      success__isnull=True,
                                      created_at__gt=JOB_CREATED_AT_CUTOFF)),
        )
        # via https://stackoverflow.com/questions/32520655/get-average-of-difference-of-datetime-fields-in-django
        result['average_time'] = jobs.filter(start_filter).filter(
            start_time__isnull=False,
            end_time__isnull=False,
            success=True
        ).aggregate(
            average_time=Avg(F('end_time') - F('start_time'))
        )['average_time']

        if not result['average_time']:
            result['average_time'] = 0
        else:
            result['average_time'] = result['average_time'].total_seconds()

        if range_param:
            result['timeline'] = cls._get_intervals(jobs, range_param) \
                                     .annotate(
                                         total=Count('id'),
                                         successful=Count('id', filter=Q(success=True)),
                                         failed=Count('id', filter=Q(success=False)),
                                         pending=Count('id', filter=Q(start_time__isnull=True, success__isnull=True)),
                                         open=Count('id', filter=Q(start_time__isnull=False, success__isnull=True)),
                                     )

        return result

    @classmethod
    def _get_object_stats(cls, objects, range_param=False, field='created_at'):
        result = {
            'total': objects.count()
        }

        if range_param:
            result['timeline'] = cls._get_intervals(objects, range_param, field)\
                                     .annotate(total=Count('id'))

        return result

    @classmethod
    def _get_intervals(cls, objects, range_param, field='last_modified'):
        range_to_trunc = {
            'day': 'hour',
            'week': 'day',
            'month': 'day',
            'year': 'month'
        }
        current_date = datetime.now(tz=timezone.utc)
        range_to_start_date = {
            'day': current_date - timedelta(days=1),
            'week': current_date - timedelta(weeks=1),
            'month': current_date - timedelta(days=30),
            'year': current_date - timedelta(days=365)
        }

        # truncate the `last_modified` field by hour, day or month depending on the `range` param
        # and annotate each object with that. This will allow us to count the number of objects
        # on each interval with a single query
        # ref https://stackoverflow.com/a/38359913/763705
        return objects.annotate(start=Trunc(field, range_to_trunc.get(range_param), output_field=DateTimeField())) \
                      .values('start') \
                      .filter(start__gte=range_to_start_date.get(range_param))


###
# Transcriptome Indices
###

class TranscriptomeIndexList(generics.ListAPIView):
    """ List all Transcriptome Indices. These are a special type of process result, necessary for processing other SRA samples. """
    serializer_class = OrganismIndexSerializer

    def get_queryset(self):
        return OrganismIndex.objects.distinct("organism", "index_type")

@method_decorator(name='get', decorator=swagger_auto_schema(manual_parameters=[
    openapi.Parameter(
        name='organism_name', in_=openapi.IN_PATH, type=openapi.TYPE_STRING,
        description="Organism name. Eg. `MUS_MUSCULUS`",
    ),
    openapi.Parameter(
        name='length', in_=openapi.IN_QUERY, type=openapi.TYPE_STRING,
        description="",
        enum=('short', 'long',),
        default='short'
    ),
]))
class TranscriptomeIndexDetail(generics.RetrieveAPIView):
    """
    Gets the S3 url associated with the organism and length, along with other metadata about
    the transcriptome index we have stored.
    """
    serializer_class = OrganismIndexSerializer

    def get_object(self):
        organism_name = self.kwargs['organism_name'].upper()
        length = self.request.query_params.get('length', 'short')

        # Get the correct organism index object, serialize it, and return it
        transcription_length = "TRANSCRIPTOME_" + length.upper()
        try:
            organism = Organism.objects.get(name=organism_name.upper())
            organism_index = OrganismIndex.objects.exclude(s3_url__exact="")\
                                .distinct("organism", "index_type")\
                                .get(organism=organism, index_type=transcription_length)
            return organism_index
        except OrganismIndex.DoesNotExist:
            raise Http404('Organism does not exists')

###
# Compendia
###

class CompendiaDetail(APIView):
    """
    A very simple modified ComputedFile endpoint which only shows Compendia results.
    """

    @swagger_auto_schema(deprecated=True)
    def get(self, request, version, format=None):

        computed_files = ComputedFile.objects.filter(is_compendia=True, is_public=True, is_qn_target=False).order_by('-created_at')

        token_id = self.request.META.get('HTTP_API_KEY', None)

        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            serializer = CompendiaWithUrlSerializer(computed_files, many=True)
        except Exception: # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            serializer = CompendiaSerializer(computed_files, many=True)

        return Response(serializer.data)


###
# QN Targets
###

class QNTargetsAvailable(generics.ListAPIView):
    """
    This is a list of all of the organisms which have available QN Targets
    """
    serializer_class = OrganismSerializer
    paginator = None

    def get_queryset(self):
        return Organism.get_objects_with_qn_targets()

@method_decorator(name='get', decorator=swagger_auto_schema(manual_parameters=[
openapi.Parameter(
    name='organism_name', in_=openapi.IN_PATH, type=openapi.TYPE_STRING,
    description="Eg `DANIO_RERIO`, `MUS_MUSCULUS`",
)], responses={404: 'QN Target not found for the given organism.'}))
class QNTargetsDetail(generics.RetrieveAPIView):
    """
    Get a detailed view of the Quantile Normalization file for an organism.
    """
    serializer_class = QNTargetSerializer

    def get_object(self):
        organism = self.kwargs['organism_name']
        organism = organism.upper().replace(" ", "_")
        try:
            organism_id = Organism.get_object_for_name(organism).id
            annotation = ComputationalResultAnnotation.objects.filter(
                data__organism_id=organism_id,
                data__is_qn=True
            ).order_by(
                '-created_at'
            ).first()
            qn_target = annotation.result.computedfile_set.first()
        except Exception:
            raise NotFound("Don't have a target for that organism!")
        if not qn_target:
            raise NotFound("Don't have a target for that organism!!")
        return qn_target

##
# Computed Files
##

class ComputedFilesList(generics.ListAPIView):
    """
    computed_files_list

    ComputedFiles are representation of files created by data-refinery processes.

    This can also be used to fetch all the compendia files we have generated with:
    ```
    GET /computed_files?is_compendia=True&is_public=True
    ```
    """
    queryset = ComputedFile.objects.all()
    serializer_class = ComputedFileListSerializer
    filter_backends = (DjangoFilterBackend, filters.OrderingFilter,)
    filterset_fields = (
        'id',
        'is_qn_target',
        'is_smashable',
        'is_qc',
        'is_compendia',
        'quant_sf_only',
        'svd_algorithm',
        'compendia_version',
        'created_at',
        'last_modified',
    )
    ordering_fields = ('id', 'created_at', 'last_modified', 'compendia_version',)
    ordering = ('-id',)

    def get_serializer_context(self):
        """
        Extra context provided to the serializer class.
        """
        serializer_context = super(ComputedFilesList, self).get_serializer_context()
        token_id = self.request.META.get('HTTP_API_KEY', None)
        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return {**serializer_context, 'token': token}
        except Exception:  # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return serializer_context

##
# Util
##

def get_client_ip(request):
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[0]
    else:
        ip = request.META.get('REMOTE_ADDR', '')
    return ip

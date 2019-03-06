from datetime import timedelta, datetime
import requests
import mailchimp3
import nomad
from typing import Dict

from django.conf import settings
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
from rest_framework import status, filters, generics
from rest_framework.exceptions import APIException
from rest_framework.exceptions import ValidationError
from rest_framework.pagination import LimitOffsetPagination
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.settings import api_settings
from rest_framework.views import APIView

from data_refinery_api.serializers import (
    ComputationalResultSerializer,
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
from data_refinery_common.utils import get_env_variable, get_active_volumes
from data_refinery_common.logging import get_and_configure_logger
from .serializers import ExperimentDocumentSerializer

logger = get_and_configure_logger(__name__)

##
# Variables
##

RUNNING_IN_CLOUD = get_env_variable("RUNNING_IN_CLOUD")
MAILCHIMP_USER = get_env_variable("MAILCHIMP_USER")
MAILCHIMP_API_KEY = get_env_variable("MAILCHIMP_API_KEY")
MAILCHIMP_LIST_ID = get_env_variable("MAILCHIMP_LIST_ID")

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
# ElasticSearch
##
from django_elasticsearch_dsl_drf.pagination import LimitOffsetPagination as ESLimitOffsetPagination

class ExperimentDocumentView(DocumentViewSet):
    """ElasticSearch powered experiment search.

    Search can be used by affixing:

        ?search=medulloblastoma
        ?id=1
        ?search=medulloblastoma&technology=microarray&has_publication=true
        ?ordering=source_first_published

    Full examples can be found in the Django-ES-DSL-DRF docs:
        https://django-elasticsearch-dsl-drf.readthedocs.io/en/0.17.1/filtering_usage_examples.html#filtering

    """

    document = ExperimentDocument
    serializer_class = ExperimentDocumentSerializer
    pagination_class = ESLimitOffsetPagination

    # Filter backends provide different functionality we want
    filter_backends = [
        FilteringFilterBackend,
        OrderingFilterBackend,
        DefaultOrderingFilterBackend,
        CompoundSearchFilterBackend,
        FacetedSearchFilterBackend
    ]

    # Primitive
    lookup_field = 'id'

    # Define search fields
    # Is this exhaustive enough?
    search_fields = (
        'title',
        'publication_title',
        'description',
        'publication_authors',
        'submitter_institution',
        'accession_code',
        'alternate_accession_code',
        'publication_doi',
        'pubmed_id',
        'sample_metadata_fields',
        'platform_names'
    )

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
        }
    }

    # Define ordering fields
    ordering_fields = {
        'id': 'id',
        'title': 'title.raw',
        'description': 'description.raw',
        'num_total_samples': 'num_total_samples',
        'num_processed_samples': 'num_processed_samples',
        'source_first_published': 'source_first_published'
    }

    # Specify default ordering
    ordering = ('-num_total_samples', 'id', 'title', 'description', '-source_first_published')

    # Facets (aka Aggregations) provide statistics about the query result set in the API response.
    # More information here: https://github.com/barseghyanartur/django-elasticsearch-dsl-drf/blob/03a3aa716db31868ca3a71340513a993741a4177/src/django_elasticsearch_dsl_drf/filter_backends/faceted_search.py#L24
    faceted_search_fields = {
        'technology': {
            'field': 'technology',
            'facet': TermsFacet,
            'enabled': True # These are enabled by default, which is more expensive but more simple.
        },
        'organism_names': {
            'field': 'organism_names',
            'facet': TermsFacet,
            'enabled': True,
            'options': {
                'size': 999999
            }
        },
        'platform_names': {
            'field': 'platform_names',
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

    Interactive filtering allows users to explore results more easily. It can be enabled using the parameter `filter_order`.
    The filter names should be sent sepparated by commas and depending on the order in which the filters are applied the
    number of samples per filter will be different.

    """

    serializer_class = ExperimentSerializer
    pagination_class = LimitOffsetPagination

    filter_backends = (NoMarkupDjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter)
    filter_class = ExperimentFilter

    # Ordering
    ordering_fields = ('total_samples_count', 'id', 'created_at', 'source_first_published', 'accession_code',)
    ordering = ('-total_samples_count',)

    def filter_samples_count(self, queryset, name, value):
        return queryset.filter(total_samples_count=value)

    # via http://www.django-rest-framework.org/api-guide/filtering/#searchfilter
    # '^' Starts-with search.
    # '=' Exact matches.
    # '@' Full-text search.
    # '$' Regex search.
    search_fields = (   'title',
                        'description',
                        'accession_code',
                        'alternate_accession_code',
                        'protocol_description',
                        'publication_title',
                        'publication_doi',
                        'publication_authors',
                        'pubmed_id',
                        'submitter_institution',
                        'experimentannotation__data',
                        # '@sample__accession_code',
                        # '@sample__platform_name',
                        # '@sample__platform_accession_code',
                        # '@sample__organism__name',
                        # '@sample__sex',
                        # '@sample__specimen_part',
                        # '@sample__disease',
                        # '@sample__compound'
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

        filter_param_names = ['organisms__name', 'technology', 'has_publication', 'platform']
        # mapping between parameter names and category names
        filter_name_map = {
            'technology': 'technology',
            'has_publication': 'publication',
            'organisms__name': 'organism',
            'platform': 'platforms'
        }

        # With interactive filtering, the filters in the last group are calculated differently, since they should stay unchanged when applied.
        # ref https://github.com/AlexsLemonade/refinebio-frontend/issues/374#issuecomment-436373470
        # This is only enabled when the parameter `filter_order` is provided (eg `filter_order=technology,platform`)
        last_filter = self.get_last_filter()
        if last_filter and last_filter in filter_param_names:
            # 1. Calculate all filters except the one in the last category
            queryset = self.search_queryset(request.query_params)
            filter_names = [f for f in filter_param_names if f != last_filter]
            response.data['filters'] = self.get_filters(queryset, filter_names)

            # 2. Calculate the filters in the last category.
            # We use a queryset built with all filters except those in the last category
            params_without_last_category = request.query_params.copy()
            params_without_last_category.pop(last_filter)
            queryset_without_last_category = self.search_queryset(params_without_last_category)
            last_category_filters = self.get_filters(queryset_without_last_category, [last_filter])
            response.data['filters'][filter_name_map[last_filter]] = last_category_filters[filter_name_map[last_filter]]
        else:
            # Otherwise calculate the filters with the search term
            response.data['filters'] = self.get_filters(self.search_queryset(), filter_param_names)

        return response

    def get_last_filter(self):
        request = self.request
        if 'filter_order' not in request.query_params:
            return False
        filter_order = request.query_params['filter_order']
        last_filter = filter_order.split(',')[-1:][0]
        # Ensure the last filter is valid and one of the applied filters
        if not last_filter or last_filter not in request.query_params:
            return False
        return last_filter

    def get_filters(self, queryset, filters_to_calculate):
        result = {
            'technology': {},
            'publication': {},
            'organism': {},
            'platforms': {}
        }

        if 'technology' in filters_to_calculate:
            # Technology
            techs = queryset.values('technology').annotate(count=Count('sample__id', distinct=True))
            for tech in techs:
                if not tech['technology'] or not tech['technology'].strip():
                    continue
                result['technology'][tech['technology']] = tech['count']

        if 'has_publication' in filters_to_calculate:
            # Publication
            pubs = queryset.values('has_publication').annotate(count=Count('sample__id', distinct=True))
            for pub in pubs:
                if pub['has_publication']:
                    result['publication']['has_publication'] = pub['count']

        if 'organisms__name' in filters_to_calculate:
            # Organisms
            organisms = queryset.values('organisms__name').annotate(count=Count('sample__id', distinct=True))
            for organism in organisms:
                # This experiment has no ExperimentOrganism-association, which is bad.
                # This information may still live on the samples though.
                if not organism['organisms__name']:
                    continue
                result['organism'][organism['organisms__name']] = organism['count']

        if 'platform' in filters_to_calculate:
            # Platforms
            platforms = queryset.values('samples__platform_name').annotate(count=Count('sample__id', distinct=True))
            for plat in platforms:
                if plat['samples__platform_name']:
                    result['platforms'][plat['samples__platform_name']] = plat['count']

        return result

    def search_queryset(self, filter_params = False):
        if filter_params:
            queryset = ExperimentFilter(filter_params, queryset=self.get_queryset()).qs
        else:
            queryset = self.get_queryset()
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

    Adding a supplying `["ALL"]` as an experiment's accession list will add all of the associated samples.

    You must also supply `email_address` with `start`, though this will never be serialized back to you.

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

        # We convert 'ALL' into the actual accession codes given
        for key in new_data['data'].keys():
            accessions = new_data['data'][key]
            if accessions == ["ALL"]:
                experiment = get_object_or_404(Experiment, accession_code=key)
                sample_codes = list(experiment.samples.filter(is_processed=True).values_list('accession_code', flat=True))
                new_data['data'][key] = sample_codes

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
            email_ccdl_ok = self.request.data.get('email_ccdl_ok', False)
            if supplied_email_address and MAILCHIMP_API_KEY and RUNNING_IN_CLOUD and email_ccdl_ok:
                try:
                    client = mailchimp3.MailChimp(mc_api=MAILCHIMP_API_KEY, mc_user=MAILCHIMP_USER)
                    data = {
                        "email_address": supplied_email_address,
                        "status": "subscribed"
                    }
                    client.lists.members.create(MAILCHIMP_LIST_ID, data)
                except mailchimp3.mailchimpclient.MailChimpError as mc_e:
                    pass # This is likely an user-already-on-list error. It's okay.
                except Exception as e:
                    # Something outside of our control has gone wrong. It's okay.
                    logger.exception("Unexpected failure trying to add user to MailChimp list.",
                            supplied_email_address=supplied_email_address,
                            mc_user=MAILCHIMP_USER
                        )

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

                if settings.RUNNING_IN_CLOUD:
                    try:
                        try:
                            remote_ip = get_client_ip(self.request)
                            city = requests.get('https://ipapi.co/' + remote_ip + '/json/', timeout=10).json()['city']
                        except Exception:
                            city = "COULD_NOT_DETERMINE"

                        new_user_text = "New user " + supplied_email_address + " from " + city + " [" + remote_ip + "] downloaded a dataset! (" + str(old_object.id) + ")"
                        webhook_url = "https://hooks.slack.com/services/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"
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

    Also accepts:
        - `dataset_id` field instead of a list of accession codes
        - `experiment_accession_code` to return the samples associated with a given experiment

    Append the pk or accession_code to the end of this URL to see a detail view.

    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        filter_dict.pop('limit', None)
        filter_dict.pop('offset', None)
        order_by = filter_dict.pop('order_by', None)
        ids = filter_dict.pop('ids', None)
        filter_by = filter_dict.pop('filter_by', None)

        if ids is not None:
            ids = [ int(x) for x in ids.split(',')]
            filter_dict['pk__in'] = ids

        experiment_accession_code = filter_dict.pop('experiment_accession_code', None)
        if experiment_accession_code:
            experiment = get_object_or_404(Experiment.objects.values('id'), accession_code=experiment_accession_code)
            filter_dict['experiments__in'] = [experiment['id']]

        accession_codes = filter_dict.pop('accession_codes', None)
        if accession_codes:
            accession_codes = accession_codes.split(',')
            filter_dict['accession_code__in'] = accession_codes

        dataset_id = filter_dict.pop('dataset_id', None)
        if dataset_id:
            dataset = get_object_or_404(Dataset, id=dataset_id)
            # Python doesn't provide a prettier way of doing this that I know about.
            filter_dict['accession_code__in'] = [item for sublist in dataset.data.values() for item in sublist]

        samples = Sample.public_objects \
            .prefetch_related('sampleannotation_set') \
            .prefetch_related('organism') \
            .prefetch_related('results') \
            .prefetch_related('results__processor') \
            .prefetch_related('results__computationalresultannotation_set') \
            .prefetch_related('results__computedfile_set') \
            .filter(**filter_dict) \
            .order_by('-is_processed') \
            .distinct()

        if order_by:
            samples = samples.order_by(order_by)

        # case insensitive search https://docs.djangoproject.com/en/2.1/ref/models/querysets/#icontains
        if filter_by:
            samples = samples.filter(   Q(title__icontains=filter_by) |
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
                                        Q(time__icontains=filter_by) |
                                        Q(sampleannotation__data__icontains=filter_by)
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
        jobs = SurveyJob.objects.filter(**filter_dict).order_by('-id')[offset:(offset + limit)]
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
        jobs = DownloaderJob.objects.filter(**filter_dict).order_by('-id')[offset: offset + limit]
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
        jobs = ProcessorJob.objects.filter(**filter_dict).order_by('-id')[offset: offset + limit]
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
        data['samples'] = self._get_object_stats(Sample.objects)
        data['experiments'] = self._get_object_stats(Experiment.objects, range_param)
        data['processed_samples'] = self._get_object_stats(Sample.processed_objects, range_param)
        data['processed_samples']['last_hour'] = self._samples_processed_last_hour()

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

        data['processed_experiments'] = self._get_object_stats(Experiment.processed_public_objects)
        data['active_volumes'] = list(get_active_volumes())

        data['dataset'] = Dataset.objects.all().aggregate(
            total=Count('id'),
            aggregated_by_experiment=Count('id', filter=Q(aggregate_by='EXPERIMENT')),
            aggregated_by_species=Count('id', filter=Q(aggregate_by='SAMPLES')),
            scale_by_none=Count('id', filter=Q(scale_by='NONE')),
            scale_by_minmax=Count('id', filter=Q(scale_by='MINMAX')),
            scale_by_standard=Count('id', filter=Q(scale_by='STANDARD')),
            scale_by_robust=Count('id', filter=Q(scale_by='ROBUST')),
        )
        # We don't save the dates when datasets are processed, but we can use
        # `last_modified`, since datasets aren't modified again after they are processed
        data['dataset']['timeline'] = self._get_intervals(
            Dataset.objects.all().filter(is_processed=True),
            range_param,
            'last_modified'
        ).annotate(
            total=Count('id'),
            total_size=Sum('size_in_bytes')
        )

        if range_param:
            data['input_data_size'] = self._get_input_data_size()
            data['output_data_size'] = self._get_output_data_size()

        try:
            nomad_stats = self._get_nomad_jobs_breakdown()
            data['nomad_running_jobs'] = nomad_stats["nomad_running_jobs"]
            data['nomad_pending_jobs'] = nomad_stats["nomad_pending_jobs"]
            data['nomad_running_jobs_by_type'] = nomad_stats["nomad_running_jobs_by_type"]
            data['nomad_pending_jobs_by_type'] = nomad_stats["nomad_pending_jobs_by_type"]
            data['nomad_running_jobs_by_volume'] = nomad_stats["nomad_running_jobs_by_volume"]
            data['nomad_pending_jobs_by_volume'] = nomad_stats["nomad_pending_jobs_by_volume"]
        except nomad.api.exceptions.BaseNomadException:
            # Nomad is not available right now, so exclude these.
            pass

        return Response(data)

    def _samples_processed_last_hour(self):
        current_date = datetime.now(tz=timezone.utc)
        start = current_date - timedelta(hours=1)
        return Sample.processed_objects.filter(created_at__range=(start, current_date)).count()

    def _aggregate_nomad_jobs_by_type(self, jobs: Dict):
        """Aggregates the pending and running job counts for each Nomad job type.

        This is accomplished by using the stats that each
        parameterized job has about its children jobs.

        `jobs` should be a response from the Nomad API's jobs endpoint.
        """
        job_types = set()
        for job in jobs:
            # Surveyor jobs don't have ids and RAM, so handle them specially.
            if job["ID"].startswith("SURVEYOR"):
                job_types.add("SURVEYOR")
            elif job["ID"] == "SMASHER" or job["ID"] == "DOWNLOADER":
                job_types.add(job["ID"])
            else:
                # Strips out the last two underscores like so:
                # SALMON_1_16384 -> SALMON
                job_type = "_".join(job["ID"].split("_")[0:-2])
                job_types.add(job_type)

        nomad_running_jobs_by_type = {}
        nomad_pending_jobs_by_type = {}
        for job_type in job_types:
            # This will count SURVEYOR_DISPATCHER jobs as SURVEYOR
            # jobs, but I think that's fine since we barely ever run
            # SURVEYOR_DISPATCHER jobs and won't need to monitor them
            # through the dashboard.
            same_jobs = [job for job in jobs if job["ID"].startswith(job_type)]

            aggregated_pending = 0
            aggregated_running = 0
            for job in same_jobs:
                children = job["JobSummary"]["Children"]
                aggregated_pending = aggregated_pending + children["Pending"]
                aggregated_running = aggregated_running + children["Running"]

            nomad_pending_jobs_by_type[job_type] = aggregated_pending
            nomad_running_jobs_by_type[job_type] = aggregated_running

        return nomad_pending_jobs_by_type, nomad_running_jobs_by_type

    def _aggregate_nomad_jobs_by_volume(self, jobs: Dict):
        """Aggregates the job counts for each EBS volume.

        This is accomplished by using the stats that each
        parameterized job has about its children jobs.

        `jobs` should be a response from the Nomad API's jobs endpoint.
        """
        volume_ids = set()
        for job in jobs:
            # These job types don't have volume indices, so we just won't count them.
            if not job["ID"].startswith("SURVEYOR") \
               and job["ID"] != "SMASHER" \
               and job["ID"] != "DOWNLOADER":
                # Strips out the volume ID like so:
                # SALMON_1_16384 -> 1
                volume_id = "_".join(job["ID"].split("_")[-2])
                volume_ids.add(volume_id)

        nomad_running_jobs_by_volume = {}
        nomad_pending_jobs_by_volume = {}
        for volume_id in volume_ids:
            if job["ID"].startswith("SURVEYOR") \
               or job["ID"] == "SMASHER" \
               or job["ID"] == "DOWNLOADER":
                continue

            def job_has_same_volume(job: Dict) -> bool:
                """Returns true if the job is on the same volume as this iteration of the loop.

                These job types don't have volume indices, so we just
                won't count them. We theoretically could try, but it
                really would be more trouble than it's worth and this
                endpoint is already going to have a hard time returning
                a response in time.
                """
                return not job["ID"].startswith("SURVEYOR") \
                    and job["ID"] != "SMASHER" \
                    and job["ID"] != "DOWNLOADER" \
                    and job["ID"].split("_")[-2] == volume_id

            jobs_with_same_volume = [job for job in jobs if job_has_same_volume(job)]

            aggregated_pending = 0
            aggregated_running = 0
            for job in jobs_with_same_volume:
                children = job["JobSummary"]["Children"]
                aggregated_pending = aggregated_pending + children["Pending"]
                aggregated_running = aggregated_running + children["Running"]

            nomad_pending_jobs_by_volume["volume_" + str(volume_id)] = aggregated_pending
            nomad_running_jobs_by_volume["volume_" + str(volume_id)] = aggregated_running

        return nomad_pending_jobs_by_volume, nomad_running_jobs_by_volume

    def _get_nomad_jobs_breakdown(self):
        nomad_host = get_env_variable("NOMAD_HOST")
        nomad_port = get_env_variable("NOMAD_PORT", "4646")
        nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

        jobs = nomad_client.jobs.get_jobs()
        parameterized_jobs = [job for job in jobs if job['ParameterizedJob']]

        nomad_pending_jobs_by_type, nomad_running_jobs_by_type = self._aggregate_nomad_jobs_by_type(parameterized_jobs)

        # To get the total jobs for running and pending, the easiest
        # AND the most efficient way is to sum up the stats we've
        # already partially summed up.
        nomad_running_jobs = 0
        for job_type, num_jobs in nomad_running_jobs_by_type.items():
            nomad_running_jobs = nomad_running_jobs + num_jobs

        nomad_pending_jobs = 0
        for job_type, num_jobs in nomad_pending_jobs_by_type.items():
            nomad_pending_jobs = nomad_pending_jobs + num_jobs

        nomad_pending_jobs_by_volume, nomad_running_jobs_by_volume = self._aggregate_nomad_jobs_by_volume(parameterized_jobs)

        return {
            "nomad_pending_jobs": nomad_pending_jobs,
            "nomad_running_jobs": nomad_running_jobs,
            "nomad_pending_jobs_by_type": nomad_pending_jobs_by_type,
            "nomad_running_jobs_by_type": nomad_running_jobs_by_type,
            "nomad_pending_jobs_by_volume": nomad_pending_jobs_by_volume,
            "nomad_running_jobs_by_volume": nomad_running_jobs_by_volume
        }

    def _get_input_data_size(self):
        total_size = OriginalFile.objects.filter(
            sample__is_processed=True # <-- SLOW
        ).aggregate(
            Sum('size_in_bytes')
        )
        return total_size['size_in_bytes__sum'] if total_size['size_in_bytes__sum'] else 0

    def _get_output_data_size(self):
        total_size = ComputedFile.public_objects.all().filter(
            s3_bucket__isnull=False,
            s3_key__isnull=True
        ).aggregate(
            Sum('size_in_bytes')
        )
        return total_size['size_in_bytes__sum'] if total_size['size_in_bytes__sum'] else 0

    def _get_job_stats(self, jobs, range_param):
        result = jobs.aggregate(
            total=Count('id'), 
            pending=Count('id', filter=Q(start_time__isnull=True)),
            completed=Count('id', filter=Q(end_time__isnull=False)),
            successful=Count('id', filter=Q(success=True)),
            open=Count('id', filter=Q(start_time__isnull=False, end_time__isnull=True, success__isnull=True)),
        )
        # via https://stackoverflow.com/questions/32520655/get-average-of-difference-of-datetime-fields-in-django
        result['average_time'] = jobs.filter(start_time__isnull=False, end_time__isnull=False, success=True).aggregate(
                average_time=Avg(F('end_time') - F('start_time')))['average_time']

        if not result['average_time']:
            result['average_time'] = 0
        else:
            result['average_time'] = result['average_time'].total_seconds()

        if range_param:
            result['timeline'] = self._get_intervals(jobs, range_param) \
                                     .annotate(
                                         total=Count('id'),
                                         completed=Count('id', filter=Q(success=True)),
                                         pending=Count('id', filter=Q(start_time__isnull=True)),
                                         failed=Count('id', filter=Q(success=False)),
                                         open=Count('id', filter=Q(success__isnull=True)),
                                     )

        return result

    def _get_object_stats(self, objects, range_param = False):
        result = {
            'total': objects.count()
        }

        if range_param:
            result['timeline'] = self._get_intervals(objects, range_param).annotate(total=Count('id'))

        return result

    def _get_intervals(self, objects, range_param, field = 'created_at'):
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

        # trucate the `created_at` field by hour, day or month depending on the `range` param
        # and annotate each object with that. This will allow us to count the number of objects
        # on each interval with a single query
        # ref https://stackoverflow.com/a/38359913/763705
        return objects.annotate(start=Trunc(field, range_to_trunc.get(range_param), output_field=DateTimeField())) \
                      .values('start') \
                      .filter(start__gte=range_to_start_date.get(range_param)) 

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
            organism = Organism.objects.get(name=params["organism"].upper())
            organism_index = (OrganismIndex.public_objects.exclude(s3_url__exact="")
                              .distinct("organism", "index_type")
                              .get(organism=organism,
                                   index_type=transcription_length))
            serializer = OrganismIndexSerializer(organism_index)
            return Response(serializer.data)
        except OrganismIndex.DoesNotExist:
            raise Http404

###
# Compendia
###

class CompendiaDetail(APIView):
    """
    A very simple modified ComputedFile endpoint which only shows Compendia results
    """

    """List all processors."""
    def get(self, request, format=None):

        computed_files = ComputedFile.objects.filter(is_compendia=True, is_public=True, is_qn_target=False).order_by('-created_at')
        serializer = CompendiaSerializer(computed_files, many=True)
        return Response(serializer.data)


###
# QN Targets
###

class QNTargetsDetail(APIView):
    """
    Quantile Normalization Targets
    """

    """List all processors."""
    def get(self, request, format=None):
        computed_files = ComputedFile.objects.filter(is_public=True, is_qn_target=True)
        serializer = QNTargetSerializer(computed_files, many=True)
        return Response(serializer.data)

##
# Computed Files
##

class ComputedFilesList(PaginatedAPIView):
    """
    """

    def get(self, request, format=None):
        filter_dict = request.query_params.dict()
        limit = max(int(filter_dict.pop('limit', 100)), 100)
        offset = int(filter_dict.pop('offset', 0))
        jobs = ComputedFile.objects.filter(**filter_dict).order_by('-id')[offset:(offset + limit)]
        serializer = ComputedFileListSerializer(jobs, many=True)
        return Response(serializer.data)

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

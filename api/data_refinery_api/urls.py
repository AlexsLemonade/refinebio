from django.conf.urls import url, include
from django.conf import settings
from django.contrib import admin
from django.urls import include, path
from rest_framework.documentation import include_docs_urls
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.urlpatterns import format_suffix_patterns
import debug_toolbar

from data_refinery_api.views import (
    SearchAndFilter,
    ExperimentList,
    ExperimentDetail,
    SampleList,
    SampleDetail,
    OrganismList,
    PlatformList,
    InstitutionList,
    SurveyJobList,
    DownloaderJobList,
    ProcessorJobList,
    ResultsList,
    ProcessorList,
    Stats,
    CreateDatasetView,
    DatasetView,
    DatasetStatsView,
    APITokenView,
    TranscriptomeIndexDetail
)

# This provides _public_ access to the /admin interface!
# Enabling this by setting DEBUG to true this will allow unauthenticated access to the admin interface.
# Very useful for debugging (since we have no User accounts), but very dangerous for prod!
class AccessUser:
    has_module_perms = has_perm = __getattr__ = lambda s, *a, **kw: True
if settings.DEBUG:
    admin.site.has_permission = lambda r: setattr(r, 'user', AccessUser()) or True

# This class provides a friendlier root API page.
class APIRoot(APIView):
    """
    Refine.bio API

    This open API provides access to all of the Refine.bio experiment, sample and result information.

    """
    def get(self, request):
        return Response({
            'experiments': reverse('experiments', request=request),
            'samples': reverse('samples', request=request),
            'organisms': reverse('organisms', request=request),
            'platforms': reverse('platforms', request=request),
            'institutions': reverse('institutions', request=request),
            'jobs': reverse('jobs', request=request),
            'stats': reverse('stats', request=request),
            'dataset': reverse('dataset_root', request=request),
            'token': reverse('token', request=request),
            'search': reverse('search', request=request)
        })

# This class provides a friendlier jobs API page.
class JobsRoot(APIView):
    """
    Jobs!
    """
    def get(self, request):
        return Response({
            'survey': reverse('survey_jobs', request=request),
            'downloader': reverse('downloader_jobs', request=request),
            'processor': reverse('processor_jobs', request=request)
        })

# This class provides a friendlier Dataset API page.
class DatasetRoot(APIView):
    """
    Use the 'create' endpoint to create a new dataset,
    then use the returned 'id' field to see and update all the fields:

    `/dataset/id-1234-1234/`

    """
    def get(self, request):
        return Response({
            'create': reverse('create_dataset', request=request)
        })

urlpatterns = [
    path('__debug__/', include(debug_toolbar.urls)),

    # Primary search and filter interface
    url(r'^search/$', SearchAndFilter.as_view(), name="search"),

    # Endpoints / Self-documentation
    url(r'^experiments/$', ExperimentList.as_view(), name="experiments"),
    url(r'^experiments/(?P<pk>.+)/$', ExperimentDetail.as_view(), name="experiments_detail"),
    url(r'^samples/$', SampleList.as_view(), name="samples"),
    url(r'^samples/(?P<pk>[0-9]+)/$', SampleDetail.as_view(), name="samples_detail"),
    url(r'^organisms/$', OrganismList.as_view(), name="organisms"),
    url(r'^platforms/$', PlatformList.as_view(), name="platforms"),
    url(r'^institutions/$', InstitutionList.as_view(), name="institutions"),
    url(r'^results/$', ResultsList.as_view(), name="results"),
    url(r'^processors/$', ProcessorList.as_view(), name="processors"),

    # Deliverables
    url(r'^dataset/$', DatasetRoot.as_view(), name="dataset_root"),
    url(r'^dataset/create/$', CreateDatasetView.as_view(), name="create_dataset"),
    url(r'^dataset/stats/(?P<id>[0-9a-f-]+)/$', DatasetStatsView.as_view(), name="dataset_stats"),
    url(r'^dataset/(?P<id>[0-9a-f-]+)/$', DatasetView.as_view(), name="dataset"),
    url(r'^token/$', APITokenView.as_view(), name="token"),
    url(r'^token/(?P<id>[0-9a-f-]+)/$', APITokenView.as_view(), name="token_id"),

    # Jobs
    url(r'^jobs/$', JobsRoot.as_view(), name="jobs"),
    url(r'^jobs/survey/$', SurveyJobList.as_view(), name="survey_jobs"),
    url(r'^jobs/downloader/$', DownloaderJobList.as_view(), name="downloader_jobs"),
    url(r'^jobs/processor/$', ProcessorJobList.as_view(), name="processor_jobs"),

    # Dashboard Driver
    url(r'^stats/$', Stats.as_view(), name="stats"),

    # Transcriptome Indices
    url(r'^transcriptome_indices', TranscriptomeIndexDetail.as_view(),
        name="transcriptome-indices"),

    # Admin
    url(r'^admin/', admin.site.urls),

    # Core API schema docs
    url(r'^docs/', include_docs_urls(title='Refine.bio API'), name="docs_schema"),

    # Root
    url(r'^$', APIRoot.as_view(), name="api_root"),
]

# This adds support explicitly typed endpoints such that appending '.json' returns that application type.
urlpatterns = format_suffix_patterns(urlpatterns)

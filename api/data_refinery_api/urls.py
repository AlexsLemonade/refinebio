from django.conf.urls import url, include
from django.conf import settings
from django.contrib import admin
from rest_framework.documentation import include_docs_urls
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.urlpatterns import format_suffix_patterns

from data_refinery_api.views import (
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
    Stats
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
            'stats': reverse('stats', request=request)
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

urlpatterns = [
    # Endpoints / Self-documentation
    url(r'^experiments/$', ExperimentList.as_view(), name="experiments"),
    url(r'^experiments/(?P<pk>[0-9]+)/$', ExperimentDetail.as_view(), name="experiments_detail"),
    url(r'^samples/$', SampleList.as_view(), name="samples"),
    url(r'^samples/(?P<pk>[0-9]+)/$', SampleDetail.as_view(), name="samples_detail"),
    url(r'^organisms/$', OrganismList.as_view(), name="organisms"),
    url(r'^platforms/$', PlatformList.as_view(), name="platforms"),
    url(r'^institutions/$', InstitutionList.as_view(), name="institutions"),

    # Jobs
    url(r'^jobs/$', JobsRoot.as_view(), name="jobs"),
    url(r'^jobs/survey/$', SurveyJobList.as_view(), name="survey_jobs"),
    url(r'^jobs/downloader/$', DownloaderJobList.as_view(), name="downloader_jobs"),
    url(r'^jobs/processor/$', ProcessorJobList.as_view(), name="processor_jobs"),

    # Dashboard Driver
    url(r'^stats/$', Stats.as_view(), name="stats"),

    # Admin
    url(r'^admin', admin.site.urls),

    # Core API schema docs
    url(r'^docs', include_docs_urls(title='Refine.bio API'), name="docs_schema"),

    # Root
    url(r'^', APIRoot.as_view(), name="api_root"),
]

# This adds support explicitly typed endpoints such that appending '.json' returns that application type.
urlpatterns = format_suffix_patterns(urlpatterns)


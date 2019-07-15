from data_refinery_api.views import ExperimentDocumentView
from rest_framework.routers import DefaultRouter
from django.conf.urls import url, include
from django.conf import settings
from django.contrib import admin
from django.urls import include, path
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.urlpatterns import format_suffix_patterns
from rest_framework import permissions
from drf_yasg.views import get_schema_view
from drf_yasg import openapi
from django.views.generic import RedirectView

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
    ComputationalResultsList,
    ProcessorList,
    Stats,
    CreateDatasetView,
    DatasetView,
    CreateApiTokenView,
    APITokenView,
    TranscriptomeIndexList,
    TranscriptomeIndexDetail,
    QNTargetsDetail,
    QNTargetsAvailable,
    CompendiaDetail,
    ComputedFilesList
)

# This provides _public_ access to the /admin interface!
# Enabling this by setting DEBUG to true this will allow unauthenticated access to the admin interface.
# Very useful for debugging (since we have no User accounts), but very dangerous for prod!
class AccessUser:
    has_module_perms = has_perm = __getattr__ = lambda s, *a, **kw: True
if settings.DEBUG:
    admin.site.has_permission = lambda r: setattr(r, 'user', AccessUser()) or True

schema_view = get_schema_view(
    openapi.Info(
        title="Refine.bio API",
        default_version='v1',
        description="""
refine.bio is a multi-organism collection of genome-wide transcriptome or gene expression data that has been obtained from publicly available repositories and uniformly processed and normalized. refine.bio allows biologists, clinicians, and machine learning researchers to search for experiments from different source repositories all in one place and build custom data sets for their questions of interest.

The swagger-ui view can be found [here](http://api.refine.bio/swagger/).

The ReDoc view can be found [here](http://api.refine.bio/).

Additional documentation can be found at [docs.refine.bio](http://docs.refine.bio/en/latest/).

### Questions/Feedback?

If you have a question or comment, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues) or send us an email at [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).
        """,
        terms_of_service="https://www.refine.bio/terms",
        contact=openapi.Contact(email="ccdl@alexslemonade.org"),
        license=openapi.License(name="BSD License"),
    ),
    public=True,
    permission_classes=(permissions.AllowAny,),
)

urlpatterns = [
    # Primary search and filter interface
    url(r'^(?P<version>v1)/search/$', ExperimentDocumentView.as_view({'get': 'list'}), name="search"),

    url(r'^(?P<version>v1)/experiments/$', ExperimentList.as_view(), name="experiments"),
    url(r'^(?P<version>v1)/experiments/(?P<accession_code>.+)/$', ExperimentDetail.as_view(), name="experiments_detail"),
    url(r'^(?P<version>v1)/samples/$', SampleList.as_view(), name="samples"),
    url(r'^(?P<version>v1)/samples/(?P<accession_code>.+)/$', SampleDetail.as_view(), name="samples_detail"),
    url(r'^(?P<version>v1)/organisms/$', OrganismList.as_view(), name="organisms"),
    url(r'^(?P<version>v1)/platforms/$', PlatformList.as_view(), name="platforms"),
    url(r'^(?P<version>v1)/institutions/$', InstitutionList.as_view(), name="institutions"),
    url(r'^(?P<version>v1)/processors/$', ProcessorList.as_view(), name="processors"),

    # Deliverables
    url(r'^(?P<version>v1)/dataset/$', CreateDatasetView.as_view(), name="create_dataset"),
    url(r'^(?P<version>v1)/dataset/(?P<id>[0-9a-f-]+)/$', DatasetView.as_view(), name="dataset"),
    url(r'^(?P<version>v1)/token/$', CreateApiTokenView.as_view(), name="token"),
    url(r'^(?P<version>v1)/token/(?P<id>[0-9a-f-]+)/$', APITokenView.as_view(), name="token_id"),

    # Jobs
    url(r'^(?P<version>v1)/jobs/survey/$', SurveyJobList.as_view(), name="survey_jobs"),
    url(r'^(?P<version>v1)/jobs/downloader/$', DownloaderJobList.as_view(), name="downloader_jobs"),
    url(r'^(?P<version>v1)/jobs/processor/$', ProcessorJobList.as_view(), name="processor_jobs"),

    # Dashboard Driver
    url(r'^(?P<version>v1)/stats/$', Stats.as_view(), name="stats"),

    # Transcriptome Indices and QN Targets
    url(r'^(?P<version>v1)/transcriptome_indices/$', TranscriptomeIndexList.as_view(), name="transcriptome-indices"),
    url(r'^(?P<version>v1)/transcriptome_indices/(?P<organism_name>.+)$', TranscriptomeIndexDetail.as_view(), name="transcriptome-indices-read"),
    url(r'^(?P<version>v1)/qn_targets/$', QNTargetsAvailable.as_view(), name="qn-targets-available"),
    url(r'^(?P<version>v1)/qn_targets/(?P<organism_name>.+)$', QNTargetsDetail.as_view(), name="qn-targets"),

    url(r'^(?P<version>v1)/compendia/$', CompendiaDetail.as_view(), name="compendia"),
    url(r'^(?P<version>v1)/computed_files/$', ComputedFilesList.as_view(), name="computed-files"),
    url(r'^(?P<version>v1)/computational_results/$', ComputationalResultsList.as_view(), name="results"),

    # Admin
    url(r'^admin/', admin.site.urls),

    # v1 api docs
    url(r'^(?P<version>v1)/swagger/$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    url(r'^(?P<version>v1)$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),

    # Redirect root urls to latest version api docs
    url(r'^swagger/$', RedirectView.as_view(url="/v1/swagger")),
    url(r'^$', RedirectView.as_view(url="/v1")),

]

# This adds support explicitly typed endpoints such that appending '.json' returns that application type.
urlpatterns = format_suffix_patterns(urlpatterns)

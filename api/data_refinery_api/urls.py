from data_refinery_api.views import ExperimentDocumentView
from rest_framework.routers import DefaultRouter
from django.conf.urls import url, include
from django.conf import settings
from django.contrib import admin
from django.urls import include, path
from rest_framework.documentation import include_docs_urls
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.reverse import reverse
from rest_framework.urlpatterns import format_suffix_patterns
from rest_framework import permissions
from drf_yasg.views import get_schema_view
from drf_yasg import openapi

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

##
# ES
##

# TODO:
# Move this to a custom router
# so we can use 'name'
# https://www.django-rest-framework.org/api-guide/routers/
router = DefaultRouter()
router.register(r'',
                ExperimentDocumentView,
                base_name='esearch',
                )
router.include_format_suffixes = False

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
    url(r'^search/$', SearchAndFilter.as_view(), name="search"),

    # Endpoints / Self-documentation
    url(r'^experiments/$', ExperimentList.as_view(), name="experiments"),
    url(r'^experiments/(?P<accession_code>.+)/$', ExperimentDetail.as_view(), name="experiments_detail"),
    url(r'^samples/$', SampleList.as_view(), name="samples"),
    url(r'^samples/(?P<accession_code>.+)/$', SampleDetail.as_view(), name="samples_detail"),
    
    url(r'^organisms/$', OrganismList.as_view(), name="organisms"),
    url(r'^platforms/$', PlatformList.as_view(), name="platforms"),
    url(r'^institutions/$', InstitutionList.as_view(), name="institutions"),
    url(r'^results/$', ResultsList.as_view(), name="results"),
    url(r'^processors/$', ProcessorList.as_view(), name="processors"),

    # Deliverables
    url(r'^dataset/create/$', CreateDatasetView.as_view(), name="create_dataset"),
    url(r'^dataset/stats/(?P<id>[0-9a-f-]+)/$', DatasetStatsView.as_view(), name="dataset_stats"),
    url(r'^dataset/(?P<id>[0-9a-f-]+)/$', DatasetView.as_view(), name="dataset"),
    url(r'^token/$', APITokenView.as_view(), name="token"),
    url(r'^token/(?P<id>[0-9a-f-]+)/$', APITokenView.as_view(), name="token_id"),

    # Jobs
    url(r'^jobs/survey/$', SurveyJobList.as_view(), name="survey_jobs"),
    url(r'^jobs/downloader/$', DownloaderJobList.as_view(), name="downloader_jobs"),
    url(r'^jobs/processor/$', ProcessorJobList.as_view(), name="processor_jobs"),

    # Dashboard Driver
    url(r'^stats/$', Stats.as_view(), name="stats"),

    # Not-publically listed APIs
    # Transcriptome Indices and QN Targets
    url(r'^transcriptome_indices', TranscriptomeIndexDetail.as_view(),
        name="transcriptome-indices"),
    url(r'^compendia', CompendiaDetail.as_view(),
        name="compendia"),
    url(r'^qn_targets_available', QNTargetsAvailable.as_view(),
        name="qn-targets-available"),
    url(r'^qn_targets', QNTargetsDetail.as_view(),
        name="qn-targets"),
    url(r'^computed_files', ComputedFilesList.as_view(),
        name="computed-files"),

    # Admin
    url(r'^admin/', admin.site.urls),

    # Core API schema docs
    url(r'^docs/', include_docs_urls(title='Refine.bio API'), name="docs_schema"),

    # ES
    url(r'^es/', ExperimentDocumentView.as_view({'get': 'list'}), name="esearch"),

    # api docs
    url(r'^swagger/$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    # Let's put redoc at the api root
    url(r'^$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc')
]

# This adds support explicitly typed endpoints such that appending '.json' returns that application type.
urlpatterns = format_suffix_patterns(urlpatterns)

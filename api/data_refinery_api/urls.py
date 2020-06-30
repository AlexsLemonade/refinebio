from django.conf import settings
from django.conf.urls import url
from django.contrib import admin
from django.urls import include, path
from django.views.generic import RedirectView
from rest_framework import permissions

from drf_yasg import openapi
from drf_yasg.views import get_schema_view

from data_refinery_api.views import (
    AboutStats,
    APITokenView,
    CompendiumResultDetailView,
    CompendiumResultListView,
    ComputationalResultDetailView,
    ComputationalResultListView,
    ComputedFileDetailView,
    ComputedFileListView,
    CreateAPITokenView,
    CreateDatasetView,
    DatasetView,
    DownloaderJobDetailView,
    DownloaderJobListView,
    ExperimentDetailView,
    ExperimentDocumentView,
    ExperimentListView,
    FailedDownloaderJobStats,
    FailedProcessorJobStats,
    InstitutionListView,
    OrganismDetailView,
    OrganismListView,
    OriginalFileDetailView,
    OriginalFileListView,
    PlatformListView,
    ProcessorDetailView,
    ProcessorJobDetailView,
    ProcessorJobListView,
    ProcessorListView,
    QNTargetsAvailable,
    QNTargetsDetailView,
    SampleDetailView,
    SampleListView,
    Stats,
    SurveyJobDetailView,
    SurveyJobListView,
    TranscriptomeIndexDetailView,
    TranscriptomeIndexListView,
)


# error handlers
def handle404error(request, exception):
    message = "The requested resource was not found on this server."
    url = "https://api.refine.bio/"

    # check to see if the 404ed request contained a version
    if not match(r"^/v[1-9]/.*", request.path):
        message = "refine.bio API resources are only available through versioned requests."

    return JsonResponse({"message": message, "docs": url, "status_code": 404,}, status=404)


def handle500error(request):
    return JsonResponse(
        {"message": "A server error occured. This has been reported.", "status_code": 500,},
        status=500,
    )


# This provides _public_ access to the /admin interface!
# Enabling this by setting DEBUG to true this will allow unauthenticated access to the admin interface.
# Very useful for debugging (since we have no User accounts), but very dangerous for prod!
class AccessUser:
    has_module_perms = has_perm = __getattr__ = lambda s, *a, **kw: True


if settings.DEBUG:
    admin.site.has_permission = lambda r: setattr(r, "user", AccessUser()) or True

schema_view = get_schema_view(
    openapi.Info(
        title="Refine.bio API",
        default_version="v1",
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
    url(
        r"^(?P<version>v1)/",
        include(
            [
                # Primary search and filter interface
                url(r"^search/$", ExperimentDocumentView.as_view({"get": "list"}), name="search"),
                url(r"^experiments/$", ExperimentListView.as_view(), name="experiments"),
                url(
                    r"^experiments/(?P<accession_code>.+)/$",
                    ExperimentDetailView.as_view(),
                    name="experiments_detail",
                ),
                url(r"^samples/$", SampleListView.as_view(), name="samples"),
                url(
                    r"^samples/(?P<accession_code>.+)/$",
                    SampleDetailView.as_view(),
                    name="samples_detail",
                ),
                url(r"^organisms/$", OrganismListView.as_view(), name="organisms"),
                url(
                    r"^organisms/(?P<name>.+)/$",
                    OrganismDetailView.as_view(),
                    name="organisms_detail",
                ),
                url(r"^platforms/$", PlatformListView.as_view(), name="platforms"),
                url(r"^institutions/$", InstitutionListView.as_view(), name="institutions"),
                url(r"^processors/$", ProcessorListView.as_view(), name="processors"),
                url(
                    r"^processors/(?P<id>[0-9a-f-]+)/$",
                    ProcessorDetailView.as_view(),
                    name="processors_detail",
                ),
                # Deliverables
                url(r"^dataset/$", CreateDatasetView.as_view(), name="create_dataset"),
                url(r"^dataset/(?P<id>[0-9a-f-]+)/$", DatasetView.as_view(), name="dataset"),
                url(r"^token/$", CreateAPITokenView.as_view(), name="token"),
                url(r"^token/(?P<id>[0-9a-f-]+)/$", APITokenView.as_view(), name="token_id"),
                # Jobs
                url(r"^jobs/survey/$", SurveyJobListView.as_view(), name="survey_jobs"),
                url(
                    r"^jobs/survey/(?P<id>[0-9a-f-]+)/$",
                    SurveyJobDetailView.as_view(),
                    name="survey_jobs_detail",
                ),
                url(r"^jobs/downloader/$", DownloaderJobListView.as_view(), name="downloader_jobs"),
                url(
                    r"^jobs/downloader/(?P<id>[0-9a-f-]+)/$",
                    DownloaderJobDetailView.as_view(),
                    name="downloader_jobs_detail",
                ),
                url(r"^jobs/processor/$", ProcessorJobListView.as_view(), name="processor_jobs"),
                url(
                    r"^jobs/processor/(?P<id>[0-9a-f-]+)/$",
                    ProcessorJobDetailView.as_view(),
                    name="processor_jobs_detail",
                ),
                # Dashboard Driver
                url(r"^stats/$", Stats.as_view(), name="stats"),
                url(
                    r"^stats/failures/downloader$",
                    FailedDownloaderJobStats.as_view(),
                    name="stats_failed_downloader",
                ),
                url(
                    r"^stats/failures/processor$",
                    FailedProcessorJobStats.as_view(),
                    name="stats_failed_processor",
                ),
                url(r"^stats-about/$", AboutStats.as_view(), name="stats_about"),
                # Transcriptome Indices
                path(
                    "transcriptome_indices/",
                    include(
                        [
                            path(
                                "",
                                TranscriptomeIndexListView.as_view(),
                                name="transcriptome_indices",
                            ),
                            path(
                                "<int:id>",
                                TranscriptomeIndexDetailView.as_view(),
                                name="transcriptome_indices_read",
                            ),
                        ]
                    ),
                ),
                # QN Targets
                url(r"^qn_targets/$", QNTargetsAvailable.as_view(), name="qn_targets_available"),
                url(
                    r"^qn_targets/(?P<organism_name>.+)$",
                    QNTargetsDetailView.as_view(),
                    name="qn_targets",
                ),
                # Computed Files
                url(r"^computed_files/$", ComputedFileListView.as_view(), name="computed_files"),
                url(
                    r"^computed_files/(?P<id>[0-9a-f-]+)/$",
                    ComputedFileDetailView.as_view(),
                    name="computed_files_detail",
                ),
                url(r"^original_files/$", OriginalFileListView.as_view(), name="original_files"),
                url(
                    r"^original_files/(?P<id>[0-9a-f-]+)/$",
                    OriginalFileDetailView.as_view(),
                    name="original_files_detail",
                ),
                url(
                    r"^computational_results/$",
                    ComputationalResultListView.as_view(),
                    name="results",
                ),
                url(
                    r"^computational_results/(?P<id>[0-9a-f-]+)/$",
                    ComputationalResultDetailView.as_view(),
                    name="results_detail",
                ),
                # Compendia
                url(r"^compendia/$", CompendiumResultListView.as_view(), name="compendium_results"),
                url(
                    r"^compendia/(?P<id>[0-9]+)/$",
                    CompendiumResultDetailView.as_view(),
                    name="compendium_result",
                ),
                # v1 api docs
                url(
                    r"^swagger/$",
                    schema_view.with_ui("swagger", cache_timeout=0),
                    name="schema_swagger_ui",
                ),
                url(r"^$", schema_view.with_ui("redoc", cache_timeout=0), name="schema_redoc"),
            ]
        ),
    ),
    # Admin
    url(r"^admin/", admin.site.urls),
    # Redirect root urls to latest version api docs
    url(r"^swagger/$", RedirectView.as_view(url="/v1/swagger")),
    url(r"^$", RedirectView.as_view(url="/v1")),
]

# handle errors
handler404 = handle404error
handler500 = handle500error

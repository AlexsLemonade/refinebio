from re import match

from django.conf import settings
from django.contrib import admin
from django.http import JsonResponse
from django.urls import include, path, re_path
from django.views.generic import RedirectView

from drf_spectacular.views import SpectacularAPIView, SpectacularRedocView, SpectacularSwaggerView

from data_refinery_api.views import (
    AboutStats,
    APITokenView,
    CompendiumResultDetailView,
    CompendiumResultListView,
    ComputationalResultDetailView,
    ComputationalResultListView,
    ComputedFileDetailView,
    ComputedFileListView,
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

    return JsonResponse(
        {
            "message": message,
            "docs": url,
            "status_code": 404,
        },
        status=404,
    )


def handle500error(request):
    return JsonResponse(
        {
            "message": "A server error occured. This has been reported.",
            "status_code": 500,
        },
        status=500,
    )


# This provides _public_ access to the /admin interface!
# Enabling this by setting DEBUG to true this will allow unauthenticated access to the admin interface.
# Very useful for debugging (since we have no User accounts), but very dangerous for prod!
class AccessUser:
    has_module_perms = has_perm = __getattr__ = lambda s, *a, **kw: True


if settings.DEBUG:
    admin.site.has_permission = lambda r: setattr(r, "user", AccessUser()) or True

urlpatterns = [
    re_path(
        r"^(?P<version>v1)/",
        include(
            [
                # Primary search and filter interface
                re_path(
                    r"^search/$", ExperimentDocumentView.as_view({"get": "list"}), name="search"
                ),
                re_path(r"^experiments/$", ExperimentListView.as_view(), name="experiments"),
                re_path(
                    r"^experiments/(?P<accession_code>.+)/$",
                    ExperimentDetailView.as_view(),
                    name="experiments_detail",
                ),
                re_path(r"^samples/$", SampleListView.as_view(), name="samples"),
                re_path(
                    r"^samples/(?P<accession_code>.+)/$",
                    SampleDetailView.as_view(),
                    name="samples_detail",
                ),
                re_path(r"^organisms/$", OrganismListView.as_view(), name="organisms"),
                re_path(
                    r"^organisms/(?P<name>.+)/$",
                    OrganismDetailView.as_view(),
                    name="organisms_detail",
                ),
                re_path(r"^platforms/$", PlatformListView.as_view(), name="platforms"),
                re_path(r"^institutions/$", InstitutionListView.as_view(), name="institutions"),
                re_path(r"^processors/$", ProcessorListView.as_view(), name="processors"),
                re_path(
                    r"^processors/(?P<id>[0-9a-f-]+)/$",
                    ProcessorDetailView.as_view(),
                    name="processors_detail",
                ),
                # Deliverables
                re_path(
                    r"^dataset/$", DatasetView.as_view({"post": "create"}), name="create_dataset"
                ),
                re_path(
                    r"^dataset/(?P<id>[0-9a-f-]+)/$",
                    DatasetView.as_view({"get": "retrieve", "put": "update"}),
                    name="dataset",
                ),
                re_path(r"^token/$", APITokenView.as_view({"post": "create"}), name="token"),
                re_path(
                    r"^token/(?P<id>[0-9a-f-]+)/$",
                    APITokenView.as_view({"get": "retrieve", "put": "update"}),
                    name="token_id",
                ),
                # Jobs
                re_path(r"^jobs/survey/$", SurveyJobListView.as_view(), name="survey_jobs"),
                re_path(
                    r"^jobs/survey/(?P<id>[0-9a-f-]+)/$",
                    SurveyJobDetailView.as_view(),
                    name="survey_jobs_detail",
                ),
                re_path(
                    r"^jobs/downloader/$", DownloaderJobListView.as_view(), name="downloader_jobs"
                ),
                re_path(
                    r"^jobs/downloader/(?P<id>[0-9a-f-]+)/$",
                    DownloaderJobDetailView.as_view(),
                    name="downloader_jobs_detail",
                ),
                re_path(
                    r"^jobs/processor/$", ProcessorJobListView.as_view(), name="processor_jobs"
                ),
                re_path(
                    r"^jobs/processor/(?P<id>[0-9a-f-]+)/$",
                    ProcessorJobDetailView.as_view(),
                    name="processor_jobs_detail",
                ),
                # Dashboard Driver
                re_path(r"^stats/$", Stats.as_view(), name="stats"),
                re_path(
                    r"^stats/failures/downloader$",
                    FailedDownloaderJobStats.as_view(),
                    name="stats_failed_downloader",
                ),
                re_path(
                    r"^stats/failures/processor$",
                    FailedProcessorJobStats.as_view(),
                    name="stats_failed_processor",
                ),
                re_path(r"^stats-about/$", AboutStats.as_view(), name="stats_about"),
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
                re_path(
                    r"^qn_targets/$", QNTargetsAvailable.as_view(), name="qn_targets_available"
                ),
                re_path(
                    r"^qn_targets/(?P<organism_name>.+)$",
                    QNTargetsDetailView.as_view(),
                    name="qn_targets",
                ),
                # Computed Files
                re_path(
                    r"^computed_files/$", ComputedFileListView.as_view(), name="computed_files"
                ),
                re_path(
                    r"^computed_files/(?P<id>[0-9a-f-]+)/$",
                    ComputedFileDetailView.as_view(),
                    name="computed_files_detail",
                ),
                re_path(
                    r"^original_files/$", OriginalFileListView.as_view(), name="original_files"
                ),
                re_path(
                    r"^original_files/(?P<id>[0-9a-f-]+)/$",
                    OriginalFileDetailView.as_view(),
                    name="original_files_detail",
                ),
                re_path(
                    r"^computational_results/$",
                    ComputationalResultListView.as_view(),
                    name="results",
                ),
                re_path(
                    r"^computational_results/(?P<id>[0-9a-f-]+)/$",
                    ComputationalResultDetailView.as_view(),
                    name="results_detail",
                ),
                # Compendia
                re_path(
                    r"^compendia/$", CompendiumResultListView.as_view(), name="compendium_results"
                ),
                re_path(
                    r"^compendia/(?P<id>[0-9]+)/$",
                    CompendiumResultDetailView.as_view(),
                    name="compendium_result",
                ),
                # v1 api docs
                re_path(r"^schema/$", SpectacularAPIView.as_view(), name="schema"),
                re_path(
                    r"^swagger/$",
                    SpectacularSwaggerView.as_view(url="/v1/schema/"),
                    name="schema_swagger_ui",
                ),
                re_path(
                    r"^$",
                    SpectacularRedocView.as_view(url="/v1/schema/"),
                    name="schema_redoc",
                ),
            ]
        ),
    ),
    # Admin
    re_path(r"^admin/", admin.site.urls),
    # Redirect root urls to latest version api docs
    re_path(r"^swagger/$", RedirectView.as_view(url="/v1/swagger")),
    re_path(r"^$", RedirectView.as_view(url="/v1")),
]

# handle errors
handler404 = handle404error
handler500 = handle500error

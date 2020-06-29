from datetime import datetime, timedelta
from re import match

from django.db.models import OuterRef
from django.http import JsonResponse, QueryDict
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework import filters, generics, serializers, status
from rest_framework.exceptions import APIException
from rest_framework.response import Response
from rest_framework.views import APIView

from django_filters.rest_framework import DjangoFilterBackend
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    APIToken,
    CompendiumResult,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    DatasetAnnotation,
    DownloaderJob,
    Experiment,
    Organism,
    OrganismIndex,
    OriginalFile,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
    SurveyJob,
)
from data_refinery_common.utils import get_active_volumes, get_nomad_jobs, get_nomad_jobs_breakdown


logger = get_and_configure_logger(__name__)


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


##
# Util
##


def get_client_ip(request):
    x_forwarded_for = request.META.get("HTTP_X_FORWARDED_FOR")
    if x_forwarded_for:
        ip = x_forwarded_for.split(",")[0]
    else:
        ip = request.META.get("REMOTE_ADDR", "")
    return ip

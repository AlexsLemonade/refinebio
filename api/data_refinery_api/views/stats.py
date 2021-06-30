##
# Contains the views Stats, FailedDownloaderJobStats, FailedProcessorJobStats, and AboutStats
##

import functools
import itertools
from datetime import datetime, timedelta

from django.conf import settings
from django.core.cache import cache
from django.db.models import Count, DateTimeField
from django.db.models.aggregates import Avg, Sum
from django.db.models.expressions import F, Q
from django.db.models.functions import Left, Trunc
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework import status
from rest_framework.pagination import LimitOffsetPagination
from rest_framework.response import Response
from rest_framework.views import APIView

import boto3
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputedFile,
    Dataset,
    DownloaderJob,
    Experiment,
    Organism,
    OriginalFile,
    ProcessorJob,
    Sample,
    SurveyJob,
)
from data_refinery_common.utils import get_env_variable

logger = get_and_configure_logger(__name__)

JOB_CREATED_AT_CUTOFF = datetime(2019, 9, 19, tzinfo=timezone.utc)

# We want to cache all stats pages for 10 minutes to reduce the load on our
# servers, because computing all of the stats is really expensive
CACHE_TIME_SECONDS = 10 * 60


def get_start_date(range_param):
    current_date = datetime.now(tz=timezone.utc)
    return {
        "day": current_date - timedelta(days=1),
        "week": current_date - timedelta(weeks=1),
        "month": current_date - timedelta(days=30),
        "year": current_date - timedelta(days=365),
    }.get(range_param)


def paginate_queryset_response(queryset, request):
    paginator = LimitOffsetPagination()
    page_items = paginator.paginate_queryset(queryset, request)
    return Response(
        data={
            "results": [x.to_dict() for x in page_items],
            "limit": paginator.limit,
            "offset": paginator.offset,
            "count": paginator.count,
        },
        status=status.HTTP_200_OK,
    )


AWS_REGION = get_env_variable(
    "AWS_REGION", "us-east-1"
)  # Default to us-east-1 if the region variable can't be found
batch = boto3.client("batch", region_name=AWS_REGION)


PENDING_STATUSES = [
    "SUBMITTED",
    "PENDING",
    "RUNNABLE",
    "STARTING",
]


def add_type(job_json: dict) -> dict:
    """Get the type for a job based on its name and add that to the job dict."""
    # A job name is user_stage_NAME_GOES_HERE_..., so we need to
    # remove the first two underscore-delimited fields, and then what
    # comes after depends on the name
    split_name = job_json["jobName"].split("_")[2:]

    # The last field of a job name is always an ID, so strip that
    split_name = split_name[:-1]

    # If the new last field is numeric, then it must be a RAM amount, so
    # strip that too. Otherwise, the last field will be part of the job
    # name for jobs without RAM amounts, so we want to keep it.
    if split_name[-1].isnumeric():
        split_name = split_name[:-1]

    job_json["type"] = "_".join(split_name)
    return job_json


def get_jobs_in_queue(batch_job_queue: str) -> list:
    """Gets all of the information for the jobs in a batch queue."""
    jobs = []
    # AWS Batch only returns one status at a time and doesn't provide a `count` or `total`.
    for status in [*PENDING_STATUSES, "RUNNING"]:
        list_jobs_dict = batch.list_jobs(jobQueue=batch_job_queue, jobStatus=status)

        jobs.extend(list_jobs_dict["jobSummaryList"])

        while "nextToken" in list_jobs_dict and list_jobs_dict["nextToken"]:
            list_jobs_dict = batch.list_jobs(
                jobQueue=batch_job_queue, jobStatus=status, nextToken=list_jobs_dict["nextToken"],
            )
            jobs.extend(list_jobs_dict["jobSummaryList"])

    return jobs


def get_batch_jobs_breakdown(force=False):
    data = {}

    if not settings.RUNNING_IN_CLOUD and not force:
        # XXX: is this the right move here? What *should* we return if we aren't
        # running in the cloud? Does this break things?
        return data

    job_queue_lists = {}
    for queue_name in settings.AWS_BATCH_QUEUE_ALL_NAMES:
        try:
            job_queue_lists[queue_name] = list(map(add_type, get_jobs_in_queue(queue_name)))
        except Exception:
            logger.exception(f"Could not get jobs for queue {queue_name}")
            raise

    all_jobs = functools.reduce(lambda acc, x: acc + x, job_queue_lists.values(), [])

    def is_pending(job):
        return job["status"] in PENDING_STATUSES

    def is_running(job):
        return job["status"] == "RUNNING"

    def count_occurrences(pred, it):
        return functools.reduce(lambda acc, x: acc + 1 if pred(x) else acc, it, 0)

    data["pending_jobs"] = count_occurrences(is_pending, all_jobs)
    data["running_jobs"] = count_occurrences(is_running, all_jobs)

    def get_job_type(job):
        return job["type"]

    # groupby must be executed on a sorted iterable
    # ref: https://docs.python.org/3.8/library/itertools.html#itertools.groupby
    # XXX: in the old code this was filtered by if the job had the key
    # "ParameterizedJob". Does Batch have an equivalent?
    sorted_jobs_by_type = sorted(all_jobs, key=get_job_type)
    # We turn this into a dict because we need to iterate it more than once
    aggregated_jobs_by_type = {
        k: list(v) for k, v in itertools.groupby(sorted_jobs_by_type, get_job_type)
    }

    data["pending_jobs_by_type"] = {
        k: count_occurrences(is_pending, v) for k, v in aggregated_jobs_by_type.items()
    }
    data["running_jobs_by_type"] = {
        k: count_occurrences(is_running, v) for k, v in aggregated_jobs_by_type.items()
    }

    # Our original data is already aggregated by queue, so we don't need to
    # re-aggregate here.

    data["pending_jobs_by_queue"] = {
        k: count_occurrences(is_pending, v) for k, v in job_queue_lists.items()
    }
    data["running_jobs_by_queue"] = {
        k: count_occurrences(is_running, v) for k, v in job_queue_lists.items()
    }

    return data


class Stats(APIView):
    """Statistics about the health of the system."""

    @swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="range",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Specify a range from which to calculate the possible options",
                enum=("day", "week", "month", "year",),
            )
        ]
    )
    # @method_decorator(cache_page(CACHE_TIME_SECONDS))
    def get(self, request, version, format=None):
        range_param = request.query_params.dict().pop("range", None)
        is_dashboard = request.query_params.dict().pop("dashboard", False)
        if is_dashboard:
            stats = Stats.calculate_dashboard_stats(range_param)
        else:
            stats = Stats.calculate_stats(range_param)

        return Response(stats)

    @classmethod
    def calculate_dashboard_stats(cls, range_param):
        """The dashboard doesn't need all of the stats, and we can
        significantly reduce the request time by only crunching the stats the
        dashboard cares about"""
        data = {}
        data["generated_on"] = timezone.now()
        data["survey_jobs"] = cls._get_job_stats(SurveyJob.objects, range_param)
        data["downloader_jobs"] = cls._get_job_stats(DownloaderJob.objects, range_param)
        data["processor_jobs"] = cls._get_job_stats(ProcessorJob.objects, range_param)
        data["experiments"] = cls._get_object_stats(Experiment.objects, range_param)

        # processed and unprocessed samples stats
        data["unprocessed_samples"] = cls._get_object_stats(
            Sample.objects.filter(is_processed=False), range_param, "last_modified"
        )
        data["processed_samples"] = cls._get_object_stats(
            Sample.processed_objects, range_param, "last_modified"
        )

        data["dataset"] = cls._get_dataset_stats(range_param)

        data.update(get_batch_jobs_breakdown())

        return data

    @classmethod
    def calculate_stats(cls, range_param):
        data = cls.calculate_dashboard_stats(range_param)

        data["processed_samples"]["last_hour"] = cls._samples_processed_last_hour()

        data["processed_samples"]["technology"] = {}
        techs = Sample.processed_objects.values("technology").annotate(count=Count("technology"))
        for tech in techs:
            if not tech["technology"] or not tech["technology"].strip():
                continue
            data["processed_samples"]["technology"][tech["technology"]] = tech["count"]

        data["processed_samples"]["organism"] = {}
        organisms = Sample.processed_objects.values("organism__name").annotate(
            count=Count("organism__name")
        )
        for organism in organisms:
            if not organism["organism__name"]:
                continue
            data["processed_samples"]["organism"][organism["organism__name"]] = organism["count"]

        data["processed_experiments"] = cls._get_object_stats(Experiment.processed_public_objects)

        if range_param:
            data["input_data_size"] = cls._get_input_data_size()
            data["output_data_size"] = cls._get_output_data_size()

        data.update(get_batch_jobs_breakdown())

        return data

    @classmethod
    def _get_dataset_stats(cls, range_param):
        """Returns stats for processed datasets"""
        result = Dataset.processed_filtered_objects.aggregate(
            total=Count("id"),
            aggregated_by_experiment=Count("id", filter=Q(aggregate_by="EXPERIMENT")),
            aggregated_by_species=Count("id", filter=Q(aggregate_by="SPECIES")),
            scale_by_none=Count("id", filter=Q(scale_by="NONE")),
            scale_by_minmax=Count("id", filter=Q(scale_by="MINMAX")),
            scale_by_standard=Count("id", filter=Q(scale_by="STANDARD")),
            scale_by_robust=Count("id", filter=Q(scale_by="ROBUST")),
        )

        if range_param:
            # We don't save the dates when datasets are processed, but we can use
            # `last_modified`, since datasets aren't modified again after they are processed
            result["timeline"] = cls._get_intervals(
                Dataset.processed_filtered_objects, range_param, "last_modified"
            ).annotate(total=Count("id"), total_size=Sum("size_in_bytes"))
        return result

    @classmethod
    def _samples_processed_last_hour(cls):
        current_date = datetime.now(tz=timezone.utc)
        start = current_date - timedelta(hours=1)
        return Sample.processed_objects.filter(last_modified__range=(start, current_date)).count()

    @classmethod
    def _get_input_data_size(cls):
        total_size = OriginalFile.objects.filter(sample__is_processed=True).aggregate(  # <-- SLOW
            Sum("size_in_bytes")
        )
        return total_size["size_in_bytes__sum"] if total_size["size_in_bytes__sum"] else 0

    @classmethod
    def _get_output_data_size(cls):
        total_size = (
            ComputedFile.public_objects.all()
            .filter(s3_bucket__isnull=False, s3_key__isnull=True)
            .aggregate(Sum("size_in_bytes"))
        )
        return total_size["size_in_bytes__sum"] if total_size["size_in_bytes__sum"] else 0

    @classmethod
    def _get_job_stats(cls, jobs, range_param):
        start_filter = Q()

        if range_param:
            start_date = get_start_date(range_param)
            start_filter = start_filter | Q(start_time__gte=start_date) | Q(start_time__isnull=True)

        result = jobs.filter(start_filter).aggregate(
            total=Count("id"),
            successful=Count("id", filter=Q(success=True)),
            failed=Count("id", filter=Q(success=False)),
            pending=Count(
                "id",
                filter=Q(
                    start_time__isnull=True,
                    success__isnull=True,
                    created_at__gt=JOB_CREATED_AT_CUTOFF,
                ),
            ),
            open=Count(
                "id",
                filter=Q(
                    start_time__isnull=False,
                    success__isnull=True,
                    created_at__gt=JOB_CREATED_AT_CUTOFF,
                ),
            ),
        )
        # via https://stackoverflow.com/questions/32520655/get-average-of-difference-of-datetime-fields-in-django
        result["average_time"] = (
            jobs.filter(start_filter)
            .filter(start_time__isnull=False, end_time__isnull=False, success=True)
            .aggregate(average_time=Avg(F("end_time") - F("start_time")))["average_time"]
        )

        if not result["average_time"]:
            result["average_time"] = 0
        else:
            result["average_time"] = result["average_time"].total_seconds()

        if range_param:
            result["timeline"] = cls._get_intervals(jobs, range_param).annotate(
                total=Count("id"),
                successful=Count("id", filter=Q(success=True)),
                failed=Count("id", filter=Q(success=False)),
                pending=Count("id", filter=Q(start_time__isnull=True, success__isnull=True)),
                open=Count("id", filter=Q(start_time__isnull=False, success__isnull=True)),
            )

        return result

    @classmethod
    def _get_object_stats(cls, objects, range_param=False, field="created_at"):
        result = {"total": objects.count()}

        if range_param:
            result["timeline"] = cls._get_intervals(objects, range_param, field).annotate(
                total=Count("id")
            )

        return result

    @classmethod
    def _get_intervals(cls, objects, range_param, field="last_modified"):
        range_to_trunc = {"day": "hour", "week": "day", "month": "day", "year": "month"}

        # truncate the parameterized field so it can be annotated by range
        # ie. each day is composed of 24 hours...
        start_trunc = Trunc(field, range_to_trunc.get(range_param), output_field=DateTimeField())

        # get the correct start time for the range
        start_range = get_start_date(range_param)

        # annotate and filter in a single query
        # ref https://stackoverflow.com/a/38359913/763705
        return objects.annotate(start=start_trunc).values("start").filter(start__gte=start_range)


class FailedDownloaderJobStats(APIView):
    @swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="range",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Specify a range from which to calculate the possible options",
                enum=("day", "week", "month", "year",),
            )
        ]
    )
    # @method_decorator(cache_page(CACHE_TIME_SECONDS))
    def get(self, request, version, format=None):
        range_param = request.query_params.dict().pop("range", "day")
        start_date = get_start_date(range_param)
        jobs = (
            DownloaderJob.objects.filter(created_at__gt=start_date)
            .annotate(reason=Left("failure_reason", 80))
            .values("reason")
            .annotate(
                job_count=Count("reason"),
                sample_count=Count(
                    "original_files__samples",
                    distinct=True,
                    filter=Q(original_files__samples__is_processed=False),
                ),
            )
            .order_by("-job_count")
        )

        return paginate_queryset_response(jobs, request)


class FailedProcessorJobStats(APIView):
    @swagger_auto_schema(
        manual_parameters=[
            openapi.Parameter(
                name="range",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_STRING,
                description="Specify a range from which to calculate the possible options",
                enum=("day", "week", "month", "year",),
            )
        ]
    )
    # @method_decorator(cache_page(CACHE_TIME_SECONDS))
    def get(self, request, version, format=None):
        range_param = request.query_params.dict().pop("range", "day")
        start_date = get_start_date(range_param)
        jobs = (
            ProcessorJob.objects.filter(created_at__gt=start_date)
            .annotate(reason=Left("failure_reason", 80))
            .values("reason")
            .annotate(
                job_count=Count("reason"),
                sample_count=Count(
                    "original_files__samples",
                    distinct=True,
                    filter=Q(original_files__samples__is_processed=False),
                ),
            )
            .order_by("-job_count")
        )

        return paginate_queryset_response(jobs, request)


class AboutStats(APIView):
    """Returns general stats for the site, used in the about page"""

    # @method_decorator(cache_page(CACHE_TIME_SECONDS))
    def get(self, request, version, format=None):
        # static values for now
        dummy = request.query_params.dict().pop("dummy", None)
        if dummy:
            # add a dummy response, calculated these on 09/25/2019
            result = {
                "samples_available": 904953 + 391022,
                "total_size_in_bytes": 832195361132962,
                "supported_organisms": 43 + 159,
                "experiments_processed": 35785 + 8661,
            }
            return Response(result)

        result = {
            "samples_available": self._get_samples_available(),
            "total_size_in_bytes": OriginalFile.objects.aggregate(total_size=Sum("size_in_bytes"))[
                "total_size"
            ],
            "supported_organisms": self._get_supported_organisms(),
            "experiments_processed": self._get_experiments_processed(),
        }
        return Response(result)

    def _get_experiments_processed(self):
        """total experiments with at least one sample processed"""
        experiments_with_sample_processed = (
            Experiment.objects.annotate(
                processed_samples_count=Count("samples", filter=Q(samples__is_processed=True)),
            )
            .filter(Q(processed_samples_count__gt=1))
            .count()
        )
        experiments_with_sample_quant = (
            ComputedFile.objects.filter(filename="quant.sf", result__samples__is_processed=False)
            .values_list("result__samples__experiments", flat=True)
            .distinct()
            .count()
        )
        return experiments_with_sample_processed + experiments_with_sample_quant

    def _get_supported_organisms(self):
        """count organisms with qn targets or that have at least one sample with quant files"""
        organisms_with_qn_targets = Organism.objects.filter(qn_target__isnull=False).count()
        organisms_without_qn_targets = (
            Organism.objects.filter(
                qn_target__isnull=True,
                sample__is_processed=False,
                sample__technology="RNA-SEQ",
                sample__results__computedfile__filename="quant.sf",
            )
            .distinct()
            .count()
        )
        return organisms_with_qn_targets + organisms_without_qn_targets

    def _get_samples_available(self):
        """count the total number of samples that are processed or that have a quant.sf file associated with them"""
        processed_samples = Sample.objects.filter(is_processed=True).count()
        unprocessed_samples_with_quant = (
            Sample.objects.filter(
                is_processed=False, technology="RNA-SEQ", results__computedfile__filename="quant.sf"
            )
            .distinct()
            .count()
        )
        return processed_samples + unprocessed_samples_with_quant

##
# Contains the views Stats, FailedDownloaderJobStats, FailedProcessorJobStats, and AboutStats
##

from datetime import datetime, timedelta

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

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

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
from data_refinery_common.utils import get_active_volumes, get_nomad_jobs_breakdown


JOB_CREATED_AT_CUTOFF = datetime(2019, 6, 5, tzinfo=timezone.utc)


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


class Stats(APIView):
    """ Statistics about the health of the system. """

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
    @method_decorator(cache_page(10 * 60))
    def get(self, request, version, format=None):
        range_param = request.query_params.dict().pop("range", None)
        cached_stats = Stats.calculate_stats(range_param)
        return Response(cached_stats)

    @classmethod
    def calculate_stats(cls, range_param):
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
        data["active_volumes"] = list(get_active_volumes())
        data["dataset"] = cls._get_dataset_stats(range_param)

        if range_param:
            data["input_data_size"] = cls._get_input_data_size()
            data["output_data_size"] = cls._get_output_data_size()

        data.update(get_nomad_jobs_breakdown())

        return data

    EMAIL_USERNAME_BLACKLIST = [
        "arielsvn",
        "cansav09",
        "d.prasad",
        "daniel.himmelstein",
        "dv.prasad991",
        "greenescientist",
        "jaclyn.n.taroni",
        "kurt.wheeler91",
        "michael.zietz",
        "miserlou",
    ]

    @classmethod
    def _get_dataset_stats(cls, range_param):
        """Returns stats for processed datasets"""
        filter_query = Q()
        for username in Stats.EMAIL_USERNAME_BLACKLIST:
            filter_query = filter_query | Q(email_address__startswith=username)
        filter_query = filter_query | Q(email_address__endswith="@alexslemonade.org")
        processed_datasets = Dataset.objects.filter(
            is_processed=True, email_address__isnull=False
        ).exclude(filter_query)
        result = processed_datasets.aggregate(
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
                processed_datasets, range_param, "last_modified"
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
    @method_decorator(cache_page(10 * 60))
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
    @method_decorator(cache_page(10 * 60))
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
    """ Returns general stats for the site, used in the about page """

    @method_decorator(cache_page(10 * 60))
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
        """ total experiments with at least one sample processed """
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
        """ count organisms with qn targets or that have at least one sample with quant files """
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
        """ count the total number of samples that are processed or that have a quant.sf file associated with them """
        processed_samples = Sample.objects.filter(is_processed=True).count()
        unprocessed_samples_with_quant = (
            Sample.objects.filter(
                is_processed=False, technology="RNA-SEQ", results__computedfile__filename="quant.sf"
            )
            .distinct()
            .count()
        )
        return processed_samples + unprocessed_samples_with_quant

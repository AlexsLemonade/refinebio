##
# Contains DatasetView
##

from collections import defaultdict

from django.core.exceptions import ValidationError
from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from rest_framework import filters, generics, mixins, serializers, viewsets
from rest_framework.exceptions import APIException

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    APIToken,
    Dataset,
    DatasetAnnotation,
    Experiment,
    Organism,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    Sample,
)

logger = get_and_configure_logger(__name__)


def get_client_ip(request):
    x_forwarded_for = request.META.get("HTTP_X_FORWARDED_FOR")
    if x_forwarded_for:
        ip = x_forwarded_for.split(",")[0]
    else:
        ip = request.META.get("REMOTE_ADDR", "")
    return ip


def experiment_has_downloadable_samples(experiment, quant_sf_only=False):
    if quant_sf_only:
        try:
            experiment = Experiment.public_objects.get(accession_code=experiment)
        except Experiment.DoesNotExist:
            return False

        samples = experiment.sample_set.filter(
            # We only want samples with a quant.sf file associated with them
            results__computedfile__filename="quant.sf",
            results__computedfile__s3_key__isnull=False,
            results__computedfile__s3_bucket__isnull=False,
        )

        if samples.count() == 0:
            return False

    else:
        try:
            experiment = Experiment.processed_public_objects.get(accession_code=experiment)
        except Experiment.DoesNotExist:
            return False

    return True


def validate_dataset(data):
    """
    Dataset validation. Each experiment should always have at least one
    sample, all samples should be downloadable, and when starting the smasher
    there should be at least one experiment.
    """
    if data.get("data") is None or type(data["data"]) != dict:
        raise serializers.ValidationError("`data` must be a dict of lists.")

    if data.get("start") and len(data["data"]) == 0:
        raise serializers.ValidationError("`data` must contain at least one experiment..")

    accessions = []
    non_downloadable_experiments = []
    for key, value in data["data"].items():
        if type(value) != list:
            raise serializers.ValidationError(
                "`data` must be a dict of lists. Problem with `" + str(key) + "`"
            )

        if len(value) < 1:
            raise serializers.ValidationError(
                "`data` must be a dict of lists, each with one or more elements. Problem with `"
                + str(key)
                + "`"
            )

        if len(value) != len(set(value)):
            raise serializers.ValidationError("Duplicate values detected in " + str(value))

        # If they want "ALL", just make sure that the experiment has at least one downloadable sample
        if value == ["ALL"]:
            if not experiment_has_downloadable_samples(
                key, quant_sf_only=data.get("quant_sf_only", False)
            ):
                non_downloadable_experiments.append(key)

        # Otherwise, we will check that all the samples they requested are downloadable
        else:
            accessions.extend(value)

    if len(non_downloadable_experiments) != 0:
        raise serializers.ValidationError(
            {
                "message": "Experiment(s) in dataset have zero downloadable samples. See `non_downloadable_experiments` for a full list",
                "non_downloadable_experiments": non_downloadable_experiments,
            }
        )

    if len(accessions) == 0:
        return

    samples = Sample.public_objects.filter(accession_code__in=accessions)
    if samples.count() != len(accessions):
        raise serializers.ValidationError(
            {
                "message": "Sample(s) in dataset do not exist on refine.bio. See `non_downloadable_samples` for a full list",
                "non_downloadable_samples": list(
                    set(accessions) - set(s.accession_code for s in samples)
                ),
            },
        )

    if data.get("quant_sf_only", False):
        samples_without_quant_sf = samples.exclude(
            # Exclude samples that have at least one uploaded quant.sf file associated with them
            results__computedfile__filename="quant.sf",
            results__computedfile__s3_key__isnull=False,
            results__computedfile__s3_bucket__isnull=False,
        )
        if samples_without_quant_sf.count() > 0:
            raise serializers.ValidationError(
                {
                    "message": "Sample(s) in dataset are missing quant.sf files. See `non_downloadable_samples` for a full list",
                    "non_downloadable_samples": [
                        s.accession_code for s in samples_without_quant_sf
                    ],
                },
            )

    else:
        unprocessed_samples = samples.exclude(is_processed=True)
        if unprocessed_samples.count() > 0:
            raise serializers.ValidationError(
                {
                    "message": "Non-downloadable sample(s) in dataset. See `non_downloadable_samples` for a full list",
                    "non_downloadable_samples": [s.accession_code for s in unprocessed_samples],
                }
            )


class DatasetDetailsExperimentSerializer(serializers.ModelSerializer):
    """ This serializer contains all of the information about an experiment needed for the download
    page
    """

    sample_metadata = serializers.ReadOnlyField(source="sample_metadata_fields")

    class Meta:
        model = Experiment
        fields = ("title", "accession_code", "organism_names", "sample_metadata", "technology")


class DatasetSerializer(serializers.ModelSerializer):
    start = serializers.NullBooleanField(required=False)
    experiments = DatasetDetailsExperimentSerializer(
        source="get_experiments", many=True, read_only=True
    )
    organism_samples = serializers.SerializerMethodField(read_only=True)
    worker_version = serializers.SerializerMethodField(read_only=True)

    def __init__(self, *args, **kwargs):
        super(DatasetSerializer, self).__init__(*args, **kwargs)

        if "context" in kwargs:
            if "request" in kwargs["context"]:
                # only include the fields `experiments` and `organism_samples` when the param `?details=true`
                # is provided. This is used on the frontend to render the downloads page
                # thanks to https://django.cowhite.com/blog/dynamically-includeexclude-fields-to-django-rest-framwork-serializers-based-on-user-requests/
                if "details" not in kwargs["context"]["request"].query_params:
                    self.fields.pop("experiments")
                    self.fields.pop("organism_samples")
                    self.fields.pop("worker_version")

            # only include the field `download_url` if a valid token is specified
            # the token lookup happens in the view.
            if "token" not in kwargs["context"]:
                self.fields.pop("download_url")

    def create(self, validated_data):
        # "start" isn't actually a field on the Dataset model, we just use it
        # on the frontend to control when the dataset gets dispatched
        if "start" in validated_data:
            validated_data.pop("start")
        return super(DatasetSerializer, self).create(validated_data)

    class Meta:
        model = Dataset
        fields = (
            "id",
            "data",
            "aggregate_by",
            "scale_by",
            "is_processing",
            "is_processed",
            "is_available",
            "has_email",
            "email_address",
            "email_ccdl_ok",
            "expires_on",
            "s3_bucket",
            "s3_key",
            "success",
            "failure_reason",
            "created_at",
            "last_modified",
            "start",
            "size_in_bytes",
            "sha1",
            "experiments",
            "organism_samples",
            "download_url",
            "quantile_normalize",
            "quant_sf_only",
            "svd_algorithm",
            "worker_version",
        )
        extra_kwargs = {
            "data": {"required": True,},
            "id": {"read_only": True,},
            "is_processing": {"read_only": True,},
            "is_processed": {"read_only": True,},
            "is_available": {"read_only": True,},
            "email_address": {"required": False, "write_only": True},
            "email_ccdl_ok": {"required": False, "write_only": True},
            "expires_on": {"read_only": True,},
            "s3_bucket": {"read_only": True,},
            "s3_key": {"read_only": True,},
            "success": {"read_only": True,},
            "failure_reason": {"read_only": True,},
            "created_at": {"read_only": True,},
            "last_modified": {"read_only": True,},
            "size_in_bytes": {"read_only": True,},
            "sha1": {"read_only": True,},
            "download_url": {"read_only": True,},
            "worker_version": {
                "read_only": True,
                "help_text": "Returns the latest version of refine.bio that was used to build this dataset.",
            },
        }

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        validate_dataset(data)
        return data

    def get_organism_samples(self, obj):
        """
        Groups the sample accession codes inside a dataset by their organisms, eg:
        { HOMO_SAPIENS: [S1, S2], DANIO: [S3] }
        Useful to avoid sending sample information on the downloads page
        """
        samples = (
            obj.get_samples()
            .prefetch_related("organism")
            .values("organism__name", "accession_code")
            .order_by("organism__name", "accession_code")
        )

        result = defaultdict(list)
        for sample in samples:
            result[sample["organism__name"]].append(sample["accession_code"])

        return result

    def get_worker_version(self, obj):
        processor_jobs = obj.processor_jobs.order_by("-created_at").values_list(
            "worker_version", flat=True
        )
        if processor_jobs:
            return processor_jobs[0]
        else:
            return None


@method_decorator(
    name="retrieve",
    decorator=swagger_auto_schema(
        operation_description="View a single Dataset.",
        manual_parameters=[
            openapi.Parameter(
                name="details",
                in_=openapi.IN_QUERY,
                type=openapi.TYPE_BOOLEAN,
                description="When set to `True`, additional fields will be included in the response with details about the experiments in the dataset. This is used mostly on the dataset page in www.refine.bio",
            )
        ],
    ),
)
@method_decorator(
    name="update",
    decorator=swagger_auto_schema(
        operation_description="""
Modify an existing Dataset.

In order to begin smashing, an activated API key must be provided in the `API-KEY` header field of the request.
To acquire and activate an API key see the documentation for the [/token](#tag/token)
endpoint.

```py
import requests
import json

params = json.dumps({
    'data': data,
    'aggregate_by': 'EXPERIMENT',
    'start': True,
    'email_address': 'refinebio@gmail.com'
})
headers = {
    'Content-Type': 'application/json',
    'API-KEY': token_id # requested from /token
}
requests.put(host + '/v1/dataset/38879729-93c8-436d-9293-b95d3f274741/', params, headers=headers)
```
"""
    ),
)
class DatasetView(
    mixins.CreateModelMixin,
    mixins.UpdateModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet,
):
    """ View and modify a single Dataset. """

    queryset = Dataset.objects.all()
    serializer_class = DatasetSerializer
    lookup_field = "id"

    def get_serializer_context(self):
        """
        Extra context provided to the serializer class.
        """
        serializer_context = super(DatasetView, self).get_serializer_context()
        token_id = self.request.META.get("HTTP_API_KEY", None)
        try:
            token = APIToken.objects.get(id=token_id, is_activated=True)
            return {**serializer_context, "token": token}
        except (APIToken.DoesNotExist, ValidationError):
            return serializer_context

    def validate_token(self):
        # Make sure we have a valid activated token.
        token_id = self.request.data.get("token_id", None)

        if not token_id:
            token_id = self.request.META.get("HTTP_API_KEY", None)

        try:
            APIToken.objects.get(id=token_id, is_activated=True)
        except (APIToken.DoesNotExist, ValidationError):
            raise serializers.ValidationError("You must provide an active API token ID")

    @staticmethod
    def convert_ALL_to_accessions(data):
        qn_organisms = Organism.get_objects_with_qn_targets()
        for key in data["data"].keys():
            accessions = data["data"][key]
            if accessions == ["ALL"]:
                experiment = get_object_or_404(Experiment, accession_code=key)

                sample_codes = list(
                    experiment.samples.filter(
                        is_processed=True, organism__in=qn_organisms
                    ).values_list("accession_code", flat=True)
                )
                data["data"][key] = sample_codes

    def validate_email_address_is_nonempty(self):
        """Check to make sure the email exists. We call this when getting ready to dispatch a dataset"""
        supplied_email_address = self.request.data.get("email_address", None)
        if supplied_email_address is None or supplied_email_address == "":
            raise serializers.ValidationError("You must provide an email address.")

    def dispatch_job(self, serializer, obj):
        processor_job = ProcessorJob()
        processor_job.pipeline_applied = "SMASHER"
        processor_job.ram_amount = 4096
        processor_job.save()

        pjda = ProcessorJobDatasetAssociation()
        pjda.processor_job = processor_job
        pjda.dataset = obj
        pjda.save()

        job_sent = False

        try:
            # Hidden method of non-dispatching for testing purposes.
            if not self.request.data.get("no_send_job", False):
                job_sent = send_job(ProcessorPipeline.SMASHER, processor_job)
            else:
                # We didn't actually send it, but we also didn't want to.
                job_sent = True
        except Exception as e:
            # Just log whatever exception happens, because the foreman wil requeue the job anyway
            logger.error(e)

        if not job_sent:
            raise APIException(
                "Unable to queue download job. Something has gone"
                " wrong and we have been notified about it."
            )

        serializer.validated_data["is_processing"] = True
        obj = serializer.save()

        # create a new dataset annotation with the information of this request
        annotation = DatasetAnnotation()
        annotation.dataset = obj
        annotation.data = {
            "start": True,
            "ip": get_client_ip(self.request),
            "user_agent": self.request.META.get("HTTP_USER_AGENT", None),
        }
        annotation.save()

    def create_or_update(self, serializer):
        """ If `start` is set, fire off the job. Otherwise just create/update the dataset"""
        data = serializer.validated_data
        DatasetView.convert_ALL_to_accessions(data)

        if data.get("start"):
            self.validate_token()
            self.validate_email_address_is_nonempty()

            obj = serializer.save()

            self.dispatch_job(serializer, obj)
        else:
            serializer.save()

    def perform_create(self, serializer):
        # Since we are creating a new dataset, there is no way it is already processed
        self.create_or_update(serializer)

    def perform_update(self, serializer):
        # Check to make sure we have not already processed the dataset
        old_object = self.get_object()
        if old_object.is_processed:
            raise serializers.ValidationError(
                "You may not update Datasets which have already been processed"
            )
        # Don't allow critical data updates to jobs that have already been submitted,
        # but do allow email address updating.
        elif old_object.is_processing:
            self.validate_email_address_is_nonempty()

            serializer.validated_data["data"] = old_object.data
            serializer.validated_data["aggregate_by"] = old_object.aggregate_by
            serializer.save()
        else:
            self.create_or_update(serializer)

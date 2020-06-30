##
# Contains CreateDatasetView and DatasetView
##

from collections import defaultdict

from django.shortcuts import get_object_or_404
from django.utils.decorators import method_decorator
from rest_framework import filters, serializers, generics

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

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


def get_client_ip(request):
    x_forwarded_for = request.META.get("HTTP_X_FORWARDED_FOR")
    if x_forwarded_for:
        ip = x_forwarded_for.split(",")[0]
    else:
        ip = request.META.get("REMOTE_ADDR", "")
    return ip


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

        try:
            if len(value) != len(set(value)):
                raise serializers.ValidationError("Duplicate values detected in " + str(value))
        except Exception as e:
            raise serializers.ValidationError("Received bad dataset data: " + str(e))

        # If they want "ALL", just make sure that the experiment has at least one downloadable sample
        if value == ["ALL"]:
            try:
                experiment = Experiment.processed_public_objects.get(accession_code=key)
            except Exception as e:
                raise serializers.ValidationError(
                    "Experiment " + key + " does not have at least one downloadable sample"
                )
        # Otherwise, we will check that all the samples they requested are downloadable
        else:
            accessions.extend(value)

    if len(accessions) > 0:
        unprocessed_samples = Sample.public_objects.filter(
            accession_code__in=accessions, is_processed=False
        )
        if unprocessed_samples.count() > 0:
            raise serializers.ValidationError(
                "Non-downloadable sample(s) '"
                + ", ".join([s.accession_code for s in unprocessed_samples])
                + "' in dataset"
            )


class CreateDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dataset
        fields = ("id", "data", "email_address", "email_ccdl_ok")

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        try:
            validate_dataset(data)
        except Exception:
            raise
        return data


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
        try:
            validate_dataset(data)
        except Exception:
            raise
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


class CreateDatasetView(generics.CreateAPIView):
    """ Creates and returns new Datasets. """

    queryset = Dataset.objects.all()
    serializer_class = CreateDatasetSerializer


@method_decorator(
    name="get",
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
    name="patch", decorator=swagger_auto_schema(auto_schema=None)
)  # partial updates not supported
@method_decorator(
    name="put",
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
class DatasetView(generics.RetrieveUpdateAPIView):
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
        except Exception:  # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            return serializer_context

    def perform_update(self, serializer):
        """ If `start` is set, fire off the job. Disables dataset data updates after that.
        """
        old_object = self.get_object()
        old_data = old_object.data
        old_aggregate = old_object.aggregate_by
        already_processing = old_object.is_processing
        new_data = serializer.validated_data

        qn_organisms = Organism.get_objects_with_qn_targets()

        # We convert 'ALL' into the actual accession codes given
        for key in new_data["data"].keys():
            accessions = new_data["data"][key]
            if accessions == ["ALL"]:
                experiment = get_object_or_404(Experiment, accession_code=key)

                sample_codes = list(
                    experiment.samples.filter(
                        is_processed=True, organism__in=qn_organisms
                    ).values_list("accession_code", flat=True)
                )
                new_data["data"][key] = sample_codes

        if old_object.is_processed:
            raise serializers.ValidationError(
                "You may not update Datasets which have already been processed"
            )
        if new_data.get("start"):

            # Make sure we have a valid activated token.
            token_id = self.request.data.get("token_id", None)

            if not token_id:
                token_id = self.request.META.get("HTTP_API_KEY", None)

            try:
                APIToken.objects.get(id=token_id, is_activated=True)
            # General APIToken.DoesNotExist or django.core.exceptions.ValidationError
            except Exception:
                raise serializers.ValidationError("You must provide an active API token ID")

            supplied_email_address = self.request.data.get("email_address", None)

            if supplied_email_address is None:
                raise serializers.ValidationError("You must provide an email address.")

            email_ccdl_ok = self.request.data.get("email_ccdl_ok", False)

            if not already_processing:
                # Create and dispatch the new job.
                processor_job = ProcessorJob()
                processor_job.pipeline_applied = "SMASHER"
                processor_job.ram_amount = 4096
                processor_job.save()

                pjda = ProcessorJobDatasetAssociation()
                pjda.processor_job = processor_job
                pjda.dataset = old_object
                pjda.save()

                job_sent = False

                obj = serializer.save()
                if obj.email_address != supplied_email_address:
                    obj.email_address = supplied_email_address
                    obj.save()
                if email_ccdl_ok:
                    obj.email_ccdl_ok = email_ccdl_ok
                    obj.save()

                try:
                    # Hidden method of non-dispatching for testing purposes.
                    if not self.request.data.get("no_send_job", False):
                        job_sent = send_job(ProcessorPipeline.SMASHER, processor_job)
                    else:
                        # We didn't actually send it, but we also didn't want to.
                        job_sent = True
                except Exception:
                    # job_sent is already false and the exception has
                    # already been logged by send_job, so nothing to
                    # do other than catch the exception.
                    pass

                if not job_sent:
                    raise APIException(
                        "Unable to queue download job. Something has gone"
                        " wrong and we have been notified about it."
                    )

                serializer.validated_data["is_processing"] = True
                obj = serializer.save()

                # create a new dataset annotation with the information of this request
                annotation = DatasetAnnotation()
                annotation.dataset = old_object
                annotation.data = {
                    "start": True,
                    "ip": get_client_ip(self.request),
                    "user_agent": self.request.META.get("HTTP_USER_AGENT", None),
                }
                annotation.save()

                return obj

        # Don't allow critical data updates to jobs that have already been submitted,
        # but do allow email address updating.
        if already_processing:
            serializer.validated_data["data"] = old_data
            serializer.validated_data["aggregate_by"] = old_aggregate
        serializer.save()

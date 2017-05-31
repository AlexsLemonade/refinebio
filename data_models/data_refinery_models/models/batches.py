import os
import urllib
from enum import Enum
from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel
from data_refinery_models.models.surveys import SurveyJob


def get_env_variable(var_name):
    try:
        return os.environ[var_name]
    except KeyError:
        error_msg = "Set the %s environment variable" % var_name
        raise ImproperlyConfigured(error_msg)


class BatchStatuses(Enum):
    """Valid values for the status field of the Batch model."""

    NEW = "NEW"
    DOWNLOADED = "DOWNLOADED"
    PROCESSED = "PROCESSED"


class Batch(TimeTrackedModel):
    """Represents a batch of data.

    The definition of a Batch is intentionally that vague. What a batch
    is will vary from source to source. It could be a single file, or
    a group of files with some kind of logical grouping such as an
    experiment.
    """

    survey_job = models.ForeignKey(SurveyJob, on_delete=models.PROTECT)
    source_type = models.CharField(max_length=256)
    size_in_bytes = models.IntegerField()
    download_url = models.CharField(max_length=4096)
    raw_format = models.CharField(max_length=256, null=True)
    processed_format = models.CharField(max_length=256, null=True)
    pipeline_required = models.CharField(max_length=256)
    platform_accession_code = models.CharField(max_length=32)
    experiment_accession_code = models.CharField(max_length=32)
    experiment_title = models.CharField(max_length=256)
    status = models.CharField(max_length=20)
    release_date = models.DateField()
    last_uploaded_date = models.DateField()
    name = models.CharField(max_length=1024)

    # This field will denote where in our system the file can be found.
    internal_location = models.CharField(max_length=256, null=True)

    # This corresponds to the organism taxonomy ID from NCBI.
    organism_id = models.IntegerField()
    # This is the organism name as it appeared in the experiment.
    organism_name = models.CharField(max_length=256)

    class Meta:
        db_table = "batches"


class BatchKeyValue(TimeTrackedModel):
    """Tracks additional fields for Batches.

    Useful for fields that would be sparsely populated if they were
    their own columns. I.e. one source may have an extra field or two
    that are worth tracking but are specific to that source.
    """

    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "batch_key_values"

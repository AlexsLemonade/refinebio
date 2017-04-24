from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel


class SurveyJob(TimeTrackedModel):
    source_type = models.CharField(max_length=256)
    success = models.NullBooleanField(null=True)

    # The start time of the query used to replicate
    replication_started_at = models.DateTimeField(null=True)

    # The end time of the query used to replicate
    replication_ended_at = models.DateTimeField(null=True)

    # The start time of the job
    start_time = models.DateTimeField(null=True)

    # The end time of the job
    end_time = models.DateTimeField(null=True)

    class Meta:
        db_table = "survey_jobs"


class SurveyJobKeyValue(TimeTrackedModel):
    """
    This table is used for tracking fields onto a SurveyJob record that
    would be sparsely populated if it was its own column.
    I.e. one source may have an extra field or two that are worth
    tracking but are specific to that source.
    """
    survey_job = models.ForeignKey(SurveyJob, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "survey_job_key_values"

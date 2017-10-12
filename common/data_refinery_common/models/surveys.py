from django.db import models
from typing import Dict
from data_refinery_common.models.base_models import TimeTrackedModel


class SurveyJob(TimeTrackedModel):
    """Records information about a Surveyor Job."""

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

    def get_properties(self) -> Dict:
        return {pair.key: pair.value for pair in self.surveyjobkeyvalue_set.all()}

    class Meta:
        db_table = "survey_jobs"


class SurveyJobKeyValue(TimeTrackedModel):
    """Tracks additional fields for SurveyJobs.

    Useful for fields that would be sparsely populated if they were
    their own columns. I.e. one source may have an extra field or two
    that are worth tracking but are specific to that source.
    """

    survey_job = models.ForeignKey(SurveyJob, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "survey_job_key_values"

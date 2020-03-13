from django.db import models


class SurveyJobKeyValue(models.Model):
    """Tracks additional fields for SurveyJobs.

    Useful for fields that would be sparsely populated if they were
    their own columns. I.e. one source may have an extra field or two
    that are worth tracking but are specific to that source.
    """

    survey_job = models.ForeignKey("SurveyJob", on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "survey_job_key_values"

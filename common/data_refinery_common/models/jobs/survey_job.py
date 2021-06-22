from typing import Dict

from django.db import models
from django.utils import timezone

from data_refinery_common.models.jobs.job_managers import (
    FailedJobsManager,
    HungJobsManager,
    LostJobsManager,
    UnqueuedJobsManager,
)


class SurveyJob(models.Model):
    """Records information about a Surveyor Job."""

    class Meta:
        db_table = "survey_jobs"

    # Managers
    objects = models.Manager()
    failed_objects = FailedJobsManager()
    hung_objects = HungJobsManager()
    lost_objects = LostJobsManager()
    unqueued_objects = UnqueuedJobsManager()

    source_type = models.CharField(max_length=256)
    success = models.BooleanField(null=True)
    no_retry = models.BooleanField(default=False)
    batch_job_id = models.CharField(max_length=256, null=True)

    # Which AWS Batch Job Queue the job was run in.
    batch_job_queue = models.CharField(max_length=100, null=True)

    ram_amount = models.IntegerField(default=1024)

    # The start time of the job
    start_time = models.DateTimeField(null=True)

    # The end time of the job
    end_time = models.DateTimeField(null=True)

    # This field represents how many times this job has been
    # retried. It starts at 0 and each time the job has to be retried
    # it will be incremented.
    num_retries = models.IntegerField(default=0)

    # This field indicates whether or not this job has been retried
    # already or not.
    retried = models.BooleanField(default=False)

    # If the job is retried, this is the id of the new job
    retried_job = models.ForeignKey("self", on_delete=models.SET_NULL, null=True)

    # This field allows jobs to specify why they failed.
    failure_reason = models.TextField(null=True)

    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(SurveyJob, self).save(*args, **kwargs)

    def get_properties(self) -> Dict:
        """ Return all associated SurveyJobKeyValues as a dict"""
        return {pair.key: pair.value for pair in self.surveyjobkeyvalue_set.all()}

    def get_accession_code(self):
        """ Return `experiment_accession_code`, the most important code."""
        try:
            kvp = self.surveyjobkeyvalue_set.get(key="experiment_accession_code")
            return kvp.value
        except Exception:
            return None

    def __str__(self):
        return "SurveyJob " + str(self.pk) + ": " + str(self.source_type)

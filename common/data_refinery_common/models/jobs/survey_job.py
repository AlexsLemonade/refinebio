from typing import Dict

from django.db import models
from django.utils import timezone

import nomad
from nomad import Nomad

from data_refinery_common.models.jobs.job_managers import (
    FailedJobsManager,
    HungJobsManager,
    LostJobsManager,
)
from data_refinery_common.utils import get_env_variable


class SurveyJob(models.Model):
    """Records information about a Surveyor Job."""

    class Meta:
        db_table = "survey_jobs"

    # Managers
    objects = models.Manager()
    failed_objects = FailedJobsManager()
    hung_objects = HungJobsManager()
    lost_objects = LostJobsManager()

    source_type = models.CharField(max_length=256)
    success = models.NullBooleanField(null=True)
    no_retry = models.BooleanField(default=False)
    nomad_job_id = models.CharField(max_length=256, null=True)

    ram_amount = models.IntegerField(default=256)

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

    def kill_nomad_job(self) -> bool:
        if not self.nomad_job_id:
            return False

        try:
            nomad_host = get_env_variable("NOMAD_HOST")
            nomad_port = get_env_variable("NOMAD_PORT", "4646")
            nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
            nomad_client.job.deregister_job(self.nomad_job_id)
        except nomad.api.exceptions.BaseNomadException:
            return False

        return True

    def __str__(self):
        return "SurveyJob " + str(self.pk) + ": " + str(self.source_type)

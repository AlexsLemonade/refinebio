from typing import Set

from django.db import models
from django.utils import timezone

import nomad
from nomad import Nomad

from data_refinery_common.models.jobs.job_managers import (
    FailedJobsManager,
    HungJobsManager,
    LostJobsManager,
)
from data_refinery_common.models.sample import Sample
from data_refinery_common.utils import get_env_variable


class ProcessorJob(models.Model):
    """Records information about running a processor."""

    class Meta:
        db_table = "processor_jobs"

        indexes = [
            models.Index(
                fields=["created_at"],
                name="processor_jobs_created_at",
                # A partial index might be better here, given our queries we don't
                # need to index the whole table. We need to update to Django 2.2
                # for this to be supported.
                # condition=Q(success=None, retried=False, no_retry=False)
                # https://github.com/AlexsLemonade/refinebio/issues/1454
            ),
        ]

    # Managers
    objects = models.Manager()
    failed_objects = FailedJobsManager()
    hung_objects = HungJobsManager()
    lost_objects = LostJobsManager()

    # This field will contain an enumerated value specifying which
    # processor pipeline was applied during the processor job.
    pipeline_applied = models.CharField(max_length=256)

    original_files = models.ManyToManyField(
        "OriginalFile", through="ProcessorJobOriginalFileAssociation"
    )
    datasets = models.ManyToManyField("DataSet", through="ProcessorJobDataSetAssociation")
    no_retry = models.BooleanField(default=False)
    abort = models.BooleanField(default=False)

    # Resources
    ram_amount = models.IntegerField(default=2048)

    # The volume index is the instance id of an AWS EC2 machine. It looks like
    # these are 19 characters, but just to be safe we'll make the max length a
    # bit higher
    volume_index = models.CharField(max_length=25, null=True)

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.BooleanField(null=True)
    nomad_job_id = models.CharField(max_length=256, null=True)

    # This field represents how many times this job has been
    # retried. It starts at 0 and each time the job has to be retried
    # it will be incremented.
    num_retries = models.IntegerField(default=0)

    # This field indicates whether or not this job has been retried
    # already or not.
    retried = models.BooleanField(default=False)

    # This point of this field is to identify which worker ran the
    # job. A few fields may actually be required or something other
    # than just an id.
    worker_id = models.CharField(max_length=256, null=True)

    # This field corresponds to the version number of the
    # data_refinery_workers project that was used to run the job.
    worker_version = models.CharField(max_length=128, null=True)

    # This field allows jobs to specify why they failed.
    failure_reason = models.TextField(null=True)

    # If the job is retried, this is the id of the new job
    retried_job = models.ForeignKey("self", on_delete=models.SET_NULL, null=True)

    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def get_samples(self) -> Set[Sample]:
        samples = set()
        for original_file in self.original_files.all():
            for sample in original_file.samples.all():
                samples.add(sample)

        return samples

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

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ProcessorJob, self).save(*args, **kwargs)

    def __str__(self):
        return "ProcessorJob " + str(self.pk) + ": " + str(self.pipeline_applied)

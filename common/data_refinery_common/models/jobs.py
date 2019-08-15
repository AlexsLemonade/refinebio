from typing import Dict, Set

from django.db import transaction
from django.db import models
from django.utils import timezone
from nomad import Nomad

from data_refinery_common.models.models import Sample, Experiment, OriginalFile
from data_refinery_common.utils import get_env_variable


class SurveyJob(models.Model):
    """Records information about a Surveyor Job."""

    class Meta:
        db_table = "survey_jobs"

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
        except:
            return None

    def kill_nomad_job(self) -> bool:
        if not self.nomad_job_id:
            return False

        try:
            nomad_host = get_env_variable("NOMAD_HOST")
            nomad_port = get_env_variable("NOMAD_PORT", "4646")
            nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
            nomad_client.job.deregister_job(self.nomad_job_id)
        except:
            return False

        return True

    def __str__(self):
        return "SurveyJob " + str(self.pk) + ": " + str(self.source_type)


class SurveyJobKeyValue(models.Model):
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


class ProcessorJob(models.Model):
    """Records information about running a processor."""

    class Meta:
        db_table = "processor_jobs"

        indexes = [
            models.Index(
                fields=['created_at'],
                name='processor_jobs_created_at',
                # A partial index might be better here, given our queries we don't
                # need to index the whole table. We need to update to Django 2.2
                # for this to be supported.
                # condition=Q(success=None, retried=False, no_retry=False)
                # https://github.com/AlexsLemonade/refinebio/issues/1454
            ),
        ]

    # This field will contain an enumerated value specifying which
    # processor pipeline was applied during the processor job.
    pipeline_applied = models.CharField(max_length=256)

    original_files = models.ManyToManyField('OriginalFile', through='ProcessorJobOriginalFileAssociation')
    datasets = models.ManyToManyField('DataSet', through='ProcessorJobDataSetAssociation')
    no_retry = models.BooleanField(default=False)

    # Resources
    ram_amount = models.IntegerField(default=2048)
    volume_index = models.CharField(max_length=3, null=True)

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)
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
    retried_job = models.ForeignKey('self', on_delete=models.SET_NULL, null=True)

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
        except:
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

class DownloaderJob(models.Model):
    """Records information about running a Downloader."""

    class Meta:
        db_table = "downloader_jobs"

        indexes = [
            models.Index(
                fields=['created_at'],
                name='downloader_jobs_created_at',
                # condition=Q(success=None, retried=False, no_retry=False)
            ),
            models.Index(fields=['worker_id']),
        ]

    # This field contains a string which corresponds to a valid
    # Downloader Task. Valid values are enumerated in:
    # data_refinery_common.job_lookup.Downloaders
    downloader_task = models.CharField(max_length=256)
    accession_code = models.CharField(max_length=256, blank=True, null=True)
    no_retry = models.BooleanField(default=False)

    original_files = models.ManyToManyField('OriginalFile', through='DownloaderJobOriginalFileAssociation')

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)
    nomad_job_id = models.CharField(max_length=256, null=True)

    # Resources
    ram_amount = models.IntegerField(default=1024)
    volume_index = models.CharField(max_length=3, null=True)

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
    retried_job = models.ForeignKey('self', on_delete=models.PROTECT, null=True)

    # If the job was recreated because the data it downloaded got
    # lost, deleted, or corrupted then this field will be true.
    # This helps prevent an infinite loop of DownloaderJob recreation.
    was_recreated = models.BooleanField(default=False)

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

        # try:
        nomad_host = get_env_variable("NOMAD_HOST")
        nomad_port = get_env_variable("NOMAD_PORT", "4646")
        nomad_client = Nomad(nomad_host, port=int(nomad_port), timeout=30)
        nomad_client.job.deregister_job(self.nomad_job_id)
        # except:
        #     return False

        return True

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(DownloaderJob, self).save(*args, **kwargs)

    def __str__(self):
        return "DownloaderJob " + str(self.pk) + ": " + str(self.downloader_task)

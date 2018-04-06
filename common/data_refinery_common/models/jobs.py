from django.db import transaction
from django.db import models
from data_refinery_common.models.base_models import TimeTrackedModel

from data_refinery_common.models.new_models import Sample, Experiment, OriginalFile


class WorkerJob(TimeTrackedModel):
    """Base model with auto created_at and updated_at fields."""

    class Meta:
        abstract = True

    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

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

class ProcessorJob(WorkerJob):
    """Records information about running a processor."""

    class Meta:
        db_table = "processor_jobs"

    # This field will contain an enumerated value specifying which
    # processor pipeline was applied during the processor job.
    pipeline_applied = models.CharField(max_length=256)

    original_files = models.ManyToManyField('OriginalFile', through='ProcessorJobOriginalFileAssociation')

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

class DownloaderJob(models.Model):
    """Records information about running a Downloader."""

    class Meta:
        db_table = "downloader_jobs"

    # This field contains a string which corresponds to a valid
    # Downloader Task. Valid values are enumerated in:
    # data_refinery_common.job_lookup.Downloaders
    downloader_task = models.CharField(max_length=256)
    accession_code = models.CharField(max_length=256, blank=True, null=True)

    original_files = models.ManyToManyField('OriginalFile', through='DownloaderJobOriginalFileAssociation')

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

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

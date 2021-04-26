from typing import Set

from django.db import models
from django.utils import timezone

from data_refinery_common.models.jobs.job_managers import (
    FailedJobsManager,
    HungJobsManager,
    LostJobsManager,
    UnqueuedJobsManager,
)
from data_refinery_common.models.sample import Sample


class DownloaderJob(models.Model):
    """Records information about running a Downloader."""

    class Meta:
        db_table = "downloader_jobs"

        indexes = [
            models.Index(
                fields=["created_at"],
                name="downloader_jobs_created_at",
                # condition=Q(success=None, retried=False, no_retry=False)
            ),
            models.Index(fields=["worker_id"]),
        ]

    # Managers
    objects = models.Manager()
    failed_objects = FailedJobsManager()
    hung_objects = HungJobsManager()
    lost_objects = LostJobsManager()
    unqueued_objects = UnqueuedJobsManager()

    # This field contains a string which corresponds to a valid
    # Downloader Task. Valid values are enumerated in:
    # data_refinery_common.job_lookup.Downloaders
    downloader_task = models.CharField(max_length=256)
    accession_code = models.CharField(max_length=256, blank=True, null=True)
    no_retry = models.BooleanField(default=False)

    original_files = models.ManyToManyField(
        "OriginalFile", through="DownloaderJobOriginalFileAssociation"
    )

    # Tracking
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)
    batch_job_id = models.CharField(max_length=256, null=True)

    # Resources
    ram_amount = models.IntegerField(default=1024)

    # The volume index is the instance id of an AWS EC2 machine. It looks like
    # these are 19 characters, but just to be safe we'll make the max length a
    # bit higher
    volume_index = models.CharField(max_length=25, null=True)

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
    retried_job = models.ForeignKey("self", on_delete=models.PROTECT, null=True)

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

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(DownloaderJob, self).save(*args, **kwargs)

    def __str__(self):
        return "DownloaderJob " + str(self.pk) + ": " + str(self.downloader_task)

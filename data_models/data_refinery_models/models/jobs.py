from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel
from data_refinery_models.models.batches import Batch


class ProcessorJob(TimeTrackedModel):
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

    # This field will contain an enumerated value specifying which processor
    # pipeline was applied during the processor job.
    pipeline_applied = models.CharField(max_length=256)

    # This field represents how many times this job has been retried. It starts
    # at 0 and each time the job has to be retried it will be incremented.
    # At some point there will probably be some code like:
    # if job.num_retries >= 3:
    #     # do a bunch of logging
    # else:
    #     # retry the job
    num_retries = models.IntegerField(default=0)

    # This point of this field is to identify which worker ran the job.
    # A few fields may actually be required or something other than just an id.
    worker_id = models.CharField(max_length=256, null=True)

    class Meta:
        db_table = "processor_jobs"


class ProcessorJobsToBatches(TimeTrackedModel):
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    processor_job = models.ForeignKey(ProcessorJob, on_delete=models.CASCADE)

    class Meta:
        db_table = "processor_jobs_to_batches"


class DownloaderJob(TimeTrackedModel):
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

    # These two fields are analagous to the fields with the same names
    # in ProcessorJob, see their descriptions for more information
    num_retries = models.IntegerField(default=0)
    worker_id = models.CharField(max_length=256, null=True)

    class Meta:
        db_table = "downloader_jobs"


class DownloaderJobsToBatches(TimeTrackedModel):
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    downloader_job = models.ForeignKey(DownloaderJob, on_delete=models.CASCADE)

    class Meta:
        db_table = "downloader_jobs_to_batches"

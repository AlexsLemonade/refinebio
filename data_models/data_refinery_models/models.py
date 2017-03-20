from django.db import models
from django.utils import timezone


class TimeTrackedModel(models.Model):
    created_at = models.DateTimeField(editable=False)
    updated_at = models.DateTimeField()

    def save(self, *args, **kwargs):
        ''' On save, update timestamps '''
        if not self.id:
            self.created_at = timezone.now()
        self.updated_at = timezone.now()
        return super(TimeTrackedModel, self).save(*args, **kwargs)

    class Meta:
        abstract = True


class SurveyJob(TimeTrackedModel):
    source_type = models.CharField(max_length=256)

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


class BatchStatuses(Enum):
    NEW = "NEW"
    DOWNLOADING = "DOWNLOADING"
    DOWNLOADED = "DOWNLOADED"
    PROCESSING = "PROCESSING"
    PROCESSED = "PROCESSED"


class Batch(TimeTrackedModel):
    survey_job = models.ForeignKey(SurveyJob)
    source_type = models.CharField(max_length=256)
    size_in_bytes = models.IntegerField()
    download_url = models.CharField(max_length=2048)
    raw_format = models.CharField(max_length=256, null=True)
    processed_format = models.CharField(max_length=256, null=True)
    pipeline_required = models.CharField(max_length=256)
    accession_code = models.CharField(max_length=32)

    # This field where denote where in our system the file can be found
    internal_location = models.CharField(max_length=256, null=True)

    # This will utilize the organism taxonomy ID from NCBI
    organism = models.IntegerField()

    status = models.CharField(max_length=10)

    class Meta:
        db_table = "batches"


# This table is used for tacking fields onto a Batch record that would
# be sparsly populated if it was its own column.
# I.e. one source may have an extra field or two that are worth tracking
# but are specific to that source.
class BatchKeyValue(TimeTrackedModel):
    batch_id = models.ForeignKey(Batch, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "batch_key_values"


class ProcessorJob(TimeTrackedModel):
    batch_id = models.ForeignKey(Batch, on_delete=models.CASCADE)
    pipeline_applied = models.CharField(max_length=256)
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)
    num_retries = models.IntegerField(default=0)

    # This point of this field is to identify what worker ran the job.
    # A few fields may actually be required or something other than just an id.
    worker_id = models.CharField(max_length=256)

    class Meta:
        db_table = "processor_jobs"


class DownloaderJob(TimeTrackedModel):
    batch_id = models.ForeignKey(Batch, on_delete=models.CASCADE)
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)
    num_retries = models.IntegerField(default=0)
    worker_id = models.CharField(max_length=256)

    class Meta:
        db_table = "downloader_jobs"

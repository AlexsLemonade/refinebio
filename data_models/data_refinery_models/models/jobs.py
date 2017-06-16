from django.db import transaction
from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel
from data_refinery_models.models.batches import Batch


class ProcessorJob(TimeTrackedModel):
    """Records information about running a processor."""

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

    @classmethod
    @transaction.atomic
    def create_job_and_relationships(cls, *args, **kwargs):
        """Inits and saves a job, its batches, and the relationships.

        Creates and saves a single processor job. For each batch
        passed in via the 'batches' keyword, saves that batch, creates
        a ProcessorJobsToBatches record, and saves that record.
        """
        batches = kwargs.pop('batches', None)
        this_job = cls(*args, **kwargs)
        if batches is None:
            raise KeyError("The 'batches' argument must be specified.")
        else:
            this_job.save()
            for batch in batches:
                batch.save()
                processor_job_to_batch = ProcessorJobsToBatches(batch=batch,
                                                                processor_job=this_job)
                processor_job_to_batch.save()

        return this_job


class ProcessorJobsToBatches(TimeTrackedModel):
    """Represents a many to many relationship.

    Maps between ProcessorJobs and Batches.
    """

    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    processor_job = models.ForeignKey(ProcessorJob, on_delete=models.CASCADE)

    class Meta:
        db_table = "processor_jobs_to_batches"


class DownloaderJob(TimeTrackedModel):
    """Records information about running a Downloader."""

    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    success = models.NullBooleanField(null=True)

    # These two fields are analagous to the fields with the same names
    # in ProcessorJob, see their descriptions for more information
    num_retries = models.IntegerField(default=0)
    worker_id = models.CharField(max_length=256, null=True)

    class Meta:
        db_table = "downloader_jobs"

    @classmethod
    @transaction.atomic
    def create_job_and_relationships(cls, *args, **kwargs):
        """Inits and saves a job, its batches, and the relationships.

        Creates and saves a single downloader job. For each batch
        passed in via the 'batches' keyword, saves that batch, creates
        a DownloaderJobsToBatches record, and saves that record.
        """
        batches = kwargs.pop('batches', None)
        this_job = cls(*args, **kwargs)
        if batches is None:
            raise KeyError("The 'batches' argument must be specified.")
        else:
            this_job.save()
            for batch in batches:
                batch.save()
                downloader_job_to_batch = DownloaderJobsToBatches(batch=batch,
                                                                  downloader_job=this_job)
                downloader_job_to_batch.save()

        return this_job


class DownloaderJobsToBatches(TimeTrackedModel):
    """Represents a many to many relationship.

    Maps between DownloaderJobs and Batches.
    """

    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    downloader_job = models.ForeignKey(DownloaderJob, on_delete=models.CASCADE)

    class Meta:
        db_table = "downloader_jobs_to_batches"

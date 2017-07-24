from django.db import transaction
from django.db import models
from data_refinery_models.models.base_models import TimeTrackedModel
from data_refinery_models.models.batches import Batch


class WorkerJob(TimeTrackedModel):
    """Base model with auto created_at and updated_at fields."""

    class Meta:
        abstract = True

    batches = models.ManyToManyField(Batch)
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

    @classmethod
    @transaction.atomic
    def create_job_and_relationships(cls, *args, **kwargs):
        """Inits and saves a job and its relationships to its batches.

        Expects keyword arguments that could be passed to the init
        method for WorkerJob, with the addition of the keyword
        argument 'batches', which must be specified as a list of Batch
        objects which have already been saved (so they have an id).
        """
        batches = kwargs.pop('batches', None)
        this_job = cls(*args, **kwargs)
        if batches is None:
            raise KeyError("The 'batches' argument must be specified.")
        else:
            this_job.save()
            this_job.batches.add(*batches)

        return this_job


class ProcessorJob(WorkerJob):
    """Records information about running a processor."""

    class Meta:
        db_table = "processor_jobs"

    # This field will contain an enumerated value specifying which
    # processor pipeline was applied during the processor job.
    pipeline_applied = models.CharField(max_length=256)


class DownloaderJob(WorkerJob):
    """Records information about running a Downloader."""

    class Meta:
        db_table = "downloader_jobs"

    # This field contains a string which corresponds to a Celery
    # task. It needs to have all the namespaces as if it were an
    # import statement because it will be used to look up the Celery
    # task on the worker machine.
    downloader_task = models.CharField(max_length=256)

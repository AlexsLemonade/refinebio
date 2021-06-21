from django.db import models


class FailedJobsManager(models.Manager):
    """
    Only returns Failed Jobs
    """

    def get_queryset(self):
        return super().get_queryset().filter(success=False, retried=False, no_retry=False)


class HungJobsManager(models.Manager):
    """
    Only returns jobs that are potentially hung
    """

    def get_queryset(self):
        return (
            super()
            .get_queryset()
            .filter(
                success=None,
                retried=False,
                no_retry=False,
                start_time__isnull=False,
                end_time=None,
                nomad_job_id__isnull=False,
            )
        )


class LostJobsManager(models.Manager):
    """
    Only returns lost jobs
    """

    def get_queryset(self):
        return (
            super()
            .get_queryset()
            .filter(success=None, retried=False, no_retry=False, start_time=None, end_time=None,)
        )

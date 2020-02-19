from django.db import models


class PublicObjectsManager(models.Manager):
    """
    Only returns objects that have is_public
    """

    def get_queryset(self):
        return super().get_queryset().filter(is_public=True)


class ProcessedObjectsManager(models.Manager):
    """
    Only returns objects that have is_processed and is_public
    """

    def get_queryset(self):
        return super().get_queryset().filter(is_processed=True, is_public=True)


class ProcessedPublicObjectsManager(models.Manager):
    """
    Only returns Experiments that are is_public and have related is_processed Samples.
    """

    def get_queryset(self):
        return super().get_queryset().filter(is_public=True, num_processed_samples__gt=0)

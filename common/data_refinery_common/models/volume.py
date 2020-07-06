from django.db import models
from django.utils import timezone


class ActiveVolumesManager(models.Manager):
    """
    Only returns objects that have is_public
    """

    def get_queryset(self):
        return super().get_queryset().filter(is_active=True)

class InactiveVolumesManager(models.Manager):
    """
    Only returns objects that have is_public
    """

    def get_queryset(self):
        return super().get_queryset().filter(is_active=False)


class Volume(models.Model):
    """ Allows us to keep track of instance volumes """

    class Meta:
        db_table = "volumes"
        base_manager_name = "objects"

    # Volume fields
    id = models.IntegerField()
    is_active = models.BooleanField(default=False)

    # Managers
    objects = models.Manager()
    active_objects = ActiveVolumesManager()
    inactive_objects = InactiveVolumesManager()

    # Common Properties
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(APIToken, self).save(*args, **kwargs)

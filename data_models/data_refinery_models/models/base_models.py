from django.db import models
from django.utils import timezone


class TimeTrackedModel(models.Model):
    """Base model with auto created_at and updated_at fields."""

    created_at = models.DateTimeField(editable=False)
    updated_at = models.DateTimeField()

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.updated_at = current_time
        return super(TimeTrackedModel, self).save(*args, **kwargs)

    class Meta:
        abstract = True

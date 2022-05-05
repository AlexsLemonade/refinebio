from django.db import models


class TimestampedModel(models.Model):
    """Base model with auto created_at and last_modified_at fields."""

    created_at = models.DateTimeField(auto_now_add=True)
    last_modified_at = models.DateTimeField(auto_now=True)

    class Meta:
        abstract = True

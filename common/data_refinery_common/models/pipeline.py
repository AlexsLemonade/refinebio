from django.contrib.postgres.fields import ArrayField
from django.db import models


class Pipeline(models.Model):
    """Pipeline that is associated with a series of ComputationalResult records."""

    name = models.CharField(max_length=255)
    steps = ArrayField(models.IntegerField(), default=list)

    class Meta:
        db_table = "pipelines"

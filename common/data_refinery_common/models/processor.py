from django.contrib.postgres.fields import JSONField
from django.db import models


class Processor(models.Model):
    """Processor associated with a certain ComputationalResult."""

    name = models.CharField(max_length=255)
    version = models.CharField(max_length=64)
    docker_image = models.CharField(max_length=255)
    environment = JSONField(default=dict)

    class Meta:
        db_table = "processors"
        unique_together = ("name", "version", "docker_image", "environment")

    def __str__(self):
        return "Processor: %s (version: %s, docker_image: %s)" % (
            self.name,
            self.version,
            self.docker_image,
        )

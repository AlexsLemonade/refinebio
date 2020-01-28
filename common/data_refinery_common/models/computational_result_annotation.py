from django.contrib.postgres.fields import JSONField
from django.db import models
from django.utils import timezone

from data_refinery_common.models.managers import PublicObjectsManager


class ComputationalResultAnnotation(models.Model):
    """ Non-standard information associated with an ComputationalResult """

    class Meta:
        db_table = "computational_result_annotations"
        base_manager_name = "public_objects"

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    result = models.ForeignKey(
        "ComputationalResult", blank=False, null=False, on_delete=models.CASCADE
    )

    # Properties
    data = JSONField(default=dict)
    is_ccdl = models.BooleanField(default=True)

    # Common Properties
    is_public = models.BooleanField(default=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ComputationalResultAnnotation, self).save(*args, **kwargs)

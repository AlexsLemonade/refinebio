from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimestampedModel
from data_refinery_common.models.managers import PublicObjectsManager


class ComputationalResultAnnotation(TimestampedModel):
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
    data = models.JSONField(default=dict)
    is_ccdl = models.BooleanField(default=True)

    # Common Properties
    is_public = models.BooleanField(default=True)

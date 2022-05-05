from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimestampedModel
from data_refinery_common.models.managers import PublicObjectsManager


class DatasetAnnotation(TimestampedModel):
    """ Semi-standard information associated with a Dataset.
    IMPORTANT: This data shouldn't not be exposed through an API. """

    class Meta:
        db_table = "dataset_annotations"
        base_manager_name = "objects"

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    dataset = models.ForeignKey("Dataset", blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = models.JSONField(default=dict)

    # Common Properties
    is_public = models.BooleanField(default=False)

from django.db import models

from data_refinery_common.models.base_models import TimestampedModel
from data_refinery_common.models.managers import PublicObjectsManager


class ExperimentAnnotation(TimestampedModel):
    """ Semi-standard information associated with an Experiment """

    class Meta:
        db_table = "experiment_annotations"
        base_manager_name = "public_objects"

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    experiment = models.ForeignKey("Experiment", blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = models.JSONField(default=dict)
    is_ccdl = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=True)

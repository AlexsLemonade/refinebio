from django.db import models
from django.utils import timezone

from data_refinery_common.models.managers import PublicObjectsManager


class ExperimentAnnotation(models.Model):
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
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ExperimentAnnotation, self).save(*args, **kwargs)

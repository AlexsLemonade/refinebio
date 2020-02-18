import uuid

from django.conf import settings
from django.db import models
from django.utils import timezone


class APIToken(models.Model):
    """ Required for starting a smash job """

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Activation
    is_activated = models.BooleanField(default=False)

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

    @property
    def terms_and_conditions(self):
        """ """
        return settings.TERMS_AND_CONDITIONS

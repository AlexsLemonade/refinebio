import uuid

from django.conf import settings
from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimestampedModel


class APIToken(TimestampedModel):
    """ Required for starting a smash job """

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Activation
    is_activated = models.BooleanField(default=False)

    @property
    def terms_and_conditions(self):
        """ """
        return settings.TERMS_AND_CONDITIONS

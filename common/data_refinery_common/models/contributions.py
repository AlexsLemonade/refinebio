from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimestampedModel


class Contribution(TimestampedModel):
    """This model represents an external contribution. It stores the source
    name and methods for a contribution so that keywords and attributes can
    track how they were imported.
    """

    # If the contribution comes from a GitHub PR, this should be github.com/<username>
    source_name = models.TextField()
    methods_url = models.URLField()

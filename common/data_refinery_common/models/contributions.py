from django.db import models
from django.utils import timezone


class Contribution(models.Model):
    """This model represents an external contribution. It stores the source
    name and methods for a contribution so that keywords and attributes can
    track how they were imported.
    """

    # If the contribution comes from a GitHub PR, this should be github.com/<username>
    source_name = models.TextField()
    methods_url = models.URLField()
    created_at = models.DateTimeField(editable=False, default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        if not self.id:
            self.created_at = timezone.now()
        return super(Contribution, self).save(*args, **kwargs)

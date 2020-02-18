from django.db import models

from data_refinery_common.models.original_file import OriginalFile
from data_refinery_common.models.sample import Sample


class OriginalFileSampleAssociation(models.Model):

    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE
    )
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "original_file_sample_associations"
        unique_together = ("original_file", "sample")

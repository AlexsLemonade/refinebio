from django.db import models

from data_refinery_common.models.computed_file import ComputedFile
from data_refinery_common.models.sample import Sample


class SampleComputedFileAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    computed_file = models.ForeignKey(
        ComputedFile, blank=False, null=False, on_delete=models.CASCADE
    )

    class Meta:
        db_table = "sample_computed_file_associations"
        unique_together = ("sample", "computed_file")

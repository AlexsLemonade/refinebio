from django.db import models

from data_refinery_common.models.computational_result import ComputationalResult
from data_refinery_common.models.sample import Sample


class SampleResultAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE
    )

    class Meta:
        db_table = "sample_result_associations"
        unique_together = ("result", "sample")

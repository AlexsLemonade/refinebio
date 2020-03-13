from django.db import models

from data_refinery_common.models.computational_result import ComputationalResult
from data_refinery_common.models.experiment import Experiment


class ExperimentResultAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE
    )

    class Meta:
        db_table = "experiment_result_associations"
        unique_together = ("result", "experiment")

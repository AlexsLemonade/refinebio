from django.db import models

from data_refinery_common.models.experiment import Experiment
from data_refinery_common.models.sample import Sample


class ExperimentSampleAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_sample_associations"
        unique_together = ("experiment", "sample")

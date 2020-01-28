from django.db import models

from data_refinery_common.models.experiment import Experiment
from data_refinery_common.models.organism import Organism


class ExperimentOrganismAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_organism_associations"
        unique_together = ("experiment", "organism")

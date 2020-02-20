from django.db import models

from data_refinery_common.models.compendium_result import CompendiumResult
from data_refinery_common.models.organism import Organism


class CompendiumResultOrganismAssociation(models.Model):

    compendium_result = models.ForeignKey(
        CompendiumResult, blank=False, null=False, on_delete=models.CASCADE
    )
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "compendium_result_organism_associations"
        unique_together = ("compendium_result", "organism")

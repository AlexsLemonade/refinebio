from django.db import models
from django.utils import timezone

from data_refinery_common.models.computed_file import ComputedFile
from data_refinery_common.models.managers import PublicObjectsManager


# Compendium Computational Result
class CompendiumResult(models.Model):
    """ Computational Result For A Compendium """

    class Meta:
        db_table = "compendium_results"
        base_manager_name = "public_objects"

    def __str__(self):
        return "CompendiumResult " + str(self.pk)

    SVD_ALGORITHM_CHOICES = (
        ("NONE", "None"),
        ("RANDOMIZED", "randomized"),
        ("ARPACK", "arpack"),
    )

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    result = models.ForeignKey(
        "ComputationalResult",
        blank=False,
        null=False,
        related_name="compendium_result",
        on_delete=models.CASCADE,
    )
    primary_organism = models.ForeignKey(
        "Organism",
        blank=False,
        null=False,
        related_name="primary_compendium_results",
        on_delete=models.CASCADE,
    )
    organisms = models.ManyToManyField(
        "Organism", related_name="compendium_results", through="CompendiumResultOrganismAssociation"
    )

    # Properties
    quant_sf_only = models.BooleanField(default=False)
    compendium_version = models.IntegerField(blank=True, null=True)
    svd_algorithm = models.CharField(
        max_length=255,
        choices=SVD_ALGORITHM_CHOICES,
        default="NONE",
        help_text="The SVD algorithm that was used to impute the compendium result.",
    )

    # Common Properties
    is_public = models.BooleanField(default=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(CompendiumResult, self).save(*args, **kwargs)

    # helper
    def get_computed_file(self):
        """ Short hand method for getting the computed file for this compendium"""
        return ComputedFile.objects.filter(result=self.result).first()

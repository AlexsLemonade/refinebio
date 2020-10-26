from django.db import models
from django.utils import timezone

from data_refinery_common.constants import CURRENT_SALMON_VERSION
from data_refinery_common.models.managers import PublicObjectsManager


class OrganismIndex(models.Model):
    """ A special type of process result, necessary for processing other SRA samples """

    class Meta:
        db_table = "organism_index"
        base_manager_name = "public_objects"

    def __str__(self):
        return (
            "OrganismIndex "
            + str(self.pk)
            + ": "
            + self.organism.name
            + " ["
            + self.index_type
            + "] - "
            + str(self.salmon_version)
        )

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    organism = models.ForeignKey("Organism", blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        "ComputationalResult", blank=False, null=False, on_delete=models.CASCADE
    )

    # ex., "TRANSCRIPTOME_LONG", "TRANSCRIPTOME_SHORT"
    index_type = models.CharField(max_length=255)

    # The name of the database (for Ensembl, should specify Main, Plants, Bacteria, etc.)
    database_name = models.CharField(max_length=255)

    # This corresponds to Ensembl's release number:
    # http://ensemblgenomes.org/info/about/release_cycle
    # Determined by hitting:
    # http://rest.ensembl.org/info/software?content-type=application/json
    release_version = models.CharField(max_length=255, default="93")

    # The name of the genome assembly used which corresponds to 'GRCh38' in:
    # ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    assembly_name = models.CharField(max_length=255, default="UNKNOWN")

    # This matters, for instance salmon 0.9.0 indexes don't work with 0.10.0
    salmon_version = models.CharField(max_length=255, default=CURRENT_SALMON_VERSION)

    # We keep the director unextracted on the shared filesystem so all
    # Salmon jobs can access it.
    absolute_directory_path = models.CharField(max_length=255, blank=True, null=True, default="")
    # Common Properties
    is_public = models.BooleanField(default=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def get_computed_file(self):
        """ Short hand method for getting the computed file for this organism index"""
        return self.result.computedfile_set.first()

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(OrganismIndex, self).save(*args, **kwargs)

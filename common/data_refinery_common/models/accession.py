from django.db import models


class Accession(models.Model):
    """Accession model."""

    class Meta:
        constraints = (
            models.UniqueConstraint(
                fields=("code", "source", "technology"), name="unique_accession"
            ),
        )
        db_table = "accessions"

    code = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    last_modified_at = models.DateTimeField(auto_now=True)
    organism = models.TextField()
    published_date = models.DateTimeField()
    sample_count = models.PositiveIntegerField(default=0)
    source = models.TextField()
    technology = models.TextField()

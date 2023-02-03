from datetime import datetime

from django.db import models
from django.utils import timezone


class GatheredAccession(models.Model):
    """Gathered accession model."""

    class Meta:
        db_table = "gathered_accessions"

    accession_code = models.TextField(unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    last_modified_at = models.DateTimeField(auto_now=True)
    organism = models.TextField()
    published_date = models.DateTimeField()
    sample_count = models.PositiveIntegerField(default=0)
    source = models.TextField()
    technology = models.TextField()

    def __eq__(self, other: object) -> bool:
        """Returns True if two objects are equal. Otherwise returns False."""
        return isinstance(other, GatheredAccession) and self.accession_code == other.accession_code

    def __hash__(self) -> int:
        """Returns accession object unique hash value."""
        return hash(self.accession_code)

    def __str__(self) -> str:
        """Returns accession default string representation."""
        return ", ".join(
            (
                self.accession_code,
                self.technology,
                self.source,
                str(self.published_date.date()),
            )
        )

    @staticmethod
    def create_from_external_entry(data, source, technology, organism=None):
        """Creates accession object from MicroArray ArrayExpress entry."""
        accession = GatheredAccession()

        accession.accession_code = (
            data.get("accession") or data.get("gse") or data.get("secondary_study_accession")
        )

        organism = data.get("organism") or data.get("scientific_name") or organism
        if organism:
            accession.organism = organism.lower()

        published_date = (
            data.get("first_public") or data.get("release_date") or data.get("submission_date")
        )
        accession.published_date = timezone.make_aware(
            datetime.strptime(published_date, "%Y-%m-%d")
        )

        accession.source = source
        accession.technology = technology

        return accession

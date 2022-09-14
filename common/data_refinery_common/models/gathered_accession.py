from datetime import datetime

from django.db import models
from django.utils import timezone


class GatheredAccession(models.Model):
    """Gathered accession model."""

    class Meta:
        db_table = "gathered_accessions"

    code = models.TextField(unique=True)
    created_at = models.DateTimeField(auto_now_add=True)
    last_modified_at = models.DateTimeField(auto_now=True)
    organism = models.TextField()
    published_date = models.DateTimeField()
    sample_count = models.PositiveIntegerField(default=0)
    source = models.TextField()
    technology = models.TextField()

    def __eq__(self, other: object) -> bool:
        """Returns True if two objects are equal. Otherwise returns False."""
        return isinstance(other, GatheredAccession) and self.code == other.code

    def __hash__(self) -> int:
        """Returns accession object unique hash value."""
        return hash(self.code)

    def __str__(self) -> str:
        """Returns accession default string representation."""
        return ", ".join((self.code, self.technology, self.source, str(self.published_date.date())))

    @staticmethod
    def create_from_ma_ae_entry(entry, organism=None):
        """Creates accession object from MicroArray ArrayExpress entry."""
        accession = GatheredAccession()
        accession.code = entry["accession"]
        accession.source = "ebi_biostudies"
        accession.technology = "microarray"

        if organism:
            accession.organism = organism
        if "release_date" in entry:
            accession.published_date = timezone.make_aware(
                datetime.strptime(entry["release_date"], "%Y-%m-%d")
            )

        return accession

    @staticmethod
    def create_from_ma_geo_entry(entry):
        """Creates accession object from MicroArray GEO meta DB entry."""
        accession = GatheredAccession()
        accession.code = entry["gse"]
        accession.source = "geo_meta_db"
        accession.technology = "microarray"

        if "organism" in entry:
            accession.organism = entry["organism"].lower()
        if "submission_date" in entry:

            accession.published_date = timezone.make_aware(
                datetime.strptime(entry["submission_date"], "%Y-%m-%d")
            )

        return accession

    @staticmethod
    def create_from_rnaseq_entry(entry):
        """Creates accession object from RNA-Seq entry."""
        accession = GatheredAccession()
        accession.code = entry["secondary_study_accession"]
        accession.source = "ebi_ena_portal"
        accession.technology = "rna-seq"

        if "scientific_name" in entry:
            accession.organism = entry["scientific_name"].lower()
        if "first_public" in entry:
            accession.published_date = timezone.make_aware(
                datetime.strptime(entry["first_public"], "%Y-%m-%d")
            )

        return accession

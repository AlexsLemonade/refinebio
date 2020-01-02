from django.conf import settings
from django.contrib.postgres.fields import ArrayField, JSONField
from django.db import models, transaction
from django.db.models import Count, DateTimeField, Prefetch
from django.db.models.expressions import F, Q
from django.utils import timezone


class SurveyedAccession(models.Model):

    accession_code = models.CharField(max_length=64, unique=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        else:
            raise AssertionError("This accession has already been surveyed!")
        return super(SurveyedAccession, self).save(*args, **kwargs)

    class Meta:
        db_table = "surveyed_accessions"


class CdfCorrectedAccession(models.Model):

    accession_code = models.CharField(max_length=64, unique=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        else:
            raise AssertionError("This accession has already been surveyed!")
        return super(CdfCorrectedAccession, self).save(*args, **kwargs)

    class Meta:
        db_table = "cdf_corrected_accessions"

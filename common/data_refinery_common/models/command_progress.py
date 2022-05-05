from django.db import models
from django.utils import timezone

from data_refinery_common.models.base_models import TimestampedModel


class SurveyedAccession(TimestampedModel):

    accession_code = models.CharField(max_length=64, unique=True)

    def save(self, *args, **kwargs):
        if self.id:
            raise AssertionError("This accession has already been surveyed!")
        return super(SurveyedAccession, self).save(*args, **kwargs)

    class Meta:
        db_table = "surveyed_accessions"


class CdfCorrectedAccession(TimestampedModel):

    accession_code = models.CharField(max_length=64, unique=True)

    def save(self, *args, **kwargs):
        if self.id:
            raise AssertionError("This accession has already been surveyed!")
        return super(CdfCorrectedAccession, self).save(*args, **kwargs)

    class Meta:
        db_table = "cdf_corrected_accessions"

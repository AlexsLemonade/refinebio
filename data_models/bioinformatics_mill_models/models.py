from django.db import models
from django.utils import timezone

# Create your models here.


class TimeTrackedModel(models.Model):
    created_at = models.DateTimeField(editable=False)
    updated_at = models.DateTimeField()

    def save(self, *args, **kwargs):
        ''' On save, update timestamps '''
        if not self.id:
            self.created_at = timezone.now()
        self.updated_at = timezone.now()
        return super(TimeTrackedModel, self).save(*args, **kwargs)

    class Meta:
        abstract = True


class Batch(TimeTrackedModel):
    source_type = models.CharField(max_length=256)
    size_in_bytes = models.IntegerField()
    download_url = models.CharField(max_length=256)
    raw_format = models.CharField(max_length=256)
    processed_format = models.CharField(max_length=256)
    processor_required = models.IntegerField()
    accession_code = models.CharField(max_length=256)

    # This field where denote where in our system the file can be found
    internal_location = models.CharField(max_length=256)

    # This will at some be a meaningful integer->organism lookup thing
    organism = models.IntegerField()

    STATUSES = (
        ("NEW", "New"),
        ("DOWNLOADED", "Downloaded"),
        ("PROCESSED", "Proccessed"),
    )
    status = models.CharField(max_length=10, choices=STATUSES)

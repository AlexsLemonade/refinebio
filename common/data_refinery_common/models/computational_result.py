from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.utils import timezone

from data_refinery_common.models.computational_result_annotation import (
    ComputationalResultAnnotation,
)
from data_refinery_common.models.computed_file import ComputedFile
from data_refinery_common.models.managers import PublicObjectsManager


class ComputationalResult(models.Model):
    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"
        base_manager_name = "public_objects"

    def __str__(self):
        processor_name_str = ""
        if self.processor:
            processor_name_str = ": " + str(self.processor.name)

        return "ComputationalResult " + str(self.pk) + processor_name_str

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    commands = ArrayField(models.TextField(), default=list)
    processor = models.ForeignKey("Processor", blank=True, null=True, on_delete=models.CASCADE)

    samples = models.ManyToManyField("Sample", through="SampleResultAssociation")

    # The Organism Index used to process the sample.
    organism_index = models.ForeignKey(
        "OrganismIndex", blank=True, null=True, on_delete=models.SET_NULL
    )

    is_ccdl = models.BooleanField(default=True)

    # Stats
    time_start = models.DateTimeField(blank=True, null=True)
    time_end = models.DateTimeField(blank=True, null=True)

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
        return super(ComputationalResult, self).save(*args, **kwargs)

    def remove_computed_files_from_s3(self):
        """ Removes all associated computed files from S3.

        Use this before deleting a computational result.
        """
        for computed_file in self.computedfile_set.all():
            computed_file.delete_s3_file()

    def get_index_length(self):
        """ Pull the index_length from one of the result annotations """
        annotations = ComputationalResultAnnotation.objects.filter(result=self)

        for annotation_json in annotations:
            if "index_length" in annotation_json.data:
                return annotation_json.data["index_length"]

        return None

    def get_quant_sf_file(self):
        return (
            ComputedFile.objects.filter(
                result=self, filename="quant.sf", s3_key__isnull=False, s3_bucket__isnull=False,
            )
            .order_by("-id")
            .first()
        )

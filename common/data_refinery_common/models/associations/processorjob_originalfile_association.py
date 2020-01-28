from django.db import models

from data_refinery_common.models.jobs.processor_job import ProcessorJob
from data_refinery_common.models.original_file import OriginalFile


class ProcessorJobOriginalFileAssociation(models.Model):

    processor_job = models.ForeignKey(
        ProcessorJob, blank=False, null=False, on_delete=models.CASCADE
    )
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE
    )

    class Meta:
        db_table = "processorjob_originalfile_associations"
        unique_together = ("processor_job", "original_file")

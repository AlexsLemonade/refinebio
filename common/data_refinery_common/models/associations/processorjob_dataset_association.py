from django.db import models

from data_refinery_common.models.dataset import Dataset
from data_refinery_common.models.jobs.processor_job import ProcessorJob


class ProcessorJobDatasetAssociation(models.Model):

    processor_job = models.ForeignKey(
        ProcessorJob, blank=False, null=False, on_delete=models.CASCADE
    )
    dataset = models.ForeignKey(Dataset, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "processorjob_dataset_associations"

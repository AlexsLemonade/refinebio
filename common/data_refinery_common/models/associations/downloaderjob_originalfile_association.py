from django.db import models

from data_refinery_common.models.jobs.downloader_job import DownloaderJob
from data_refinery_common.models.original_file import OriginalFile


class DownloaderJobOriginalFileAssociation(models.Model):

    downloader_job = models.ForeignKey(
        DownloaderJob, blank=False, null=False, on_delete=models.CASCADE
    )
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE
    )

    class Meta:
        db_table = "downloaderjob_originalfile_associations"
        unique_together = ("downloader_job", "original_file")

import os

from django.conf import settings
from django.db import models
from django.utils import timezone

from data_refinery_common.constants import CURRENT_SALMON_VERSION, SYSTEM_VERSION
from data_refinery_common.enums import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.managers import PublicObjectsManager
from data_refinery_common.utils import FileUtils, calculate_file_size, calculate_md5, calculate_sha1

logger = get_and_configure_logger(__name__)


class OriginalFile(models.Model):
    """ A representation of a file from an external source """

    class Meta:
        db_table = "original_files"

        indexes = [
            models.Index(fields=["filename"]),
            models.Index(fields=["source_filename"]),
        ]

    def __str__(self):
        return "OriginalFile: " + self.get_display_name()

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # File Properties
    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField(blank=True, null=True)
    sha1 = models.CharField(max_length=64)
    md5 = models.CharField(null=True, max_length=32)
    expected_md5 = models.CharField(null=True, max_length=32)
    expected_size_in_bytes = models.BigIntegerField(blank=True, null=True)

    # AWS
    s3_bucket = models.CharField(max_length=255, blank=True, null=True)
    s3_key = models.CharField(max_length=255, blank=True, null=True)

    # Relations
    samples = models.ManyToManyField("Sample", through="OriginalFileSampleAssociation")
    processor_jobs = models.ManyToManyField(
        "data_refinery_common.ProcessorJob", through="ProcessorJobOriginalFileAssociation"
    )
    downloader_jobs = models.ManyToManyField(
        "data_refinery_common.DownloaderJob", through="DownloaderJobOriginalFileAssociation"
    )

    # Historical Properties
    source_url = models.TextField()
    is_archive = models.BooleanField(default=True)
    source_filename = models.CharField(max_length=255, blank=False)

    # Scientific Properties
    has_raw = models.BooleanField(default=True)  # Did this sample have a raw data source?

    # Crunch Properties
    is_downloaded = models.BooleanField(default=False)

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
        return super(OriginalFile, self).save(*args, **kwargs)

    def set_downloaded(self, absolute_file_path, filename=None):
        """ Marks the file as downloaded, if `filename` is not provided it will
        be parsed from the `absolute_file_path` """
        self.is_downloaded = True
        self.is_archive = FileUtils.is_archive(absolute_file_path)
        self.absolute_file_path = absolute_file_path
        self.filename = filename if filename else os.path.basename(absolute_file_path)
        self.calculate_size()
        self.calculate_sha1()
        self.calculate_md5()
        self.save()

    def calculate_sha1(self) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        self.sha1 = calculate_sha1(self.absolute_file_path)
        return self.sha1

    def calculate_md5(self) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        self.md5 = calculate_md5(self.absolute_file_path)
        return self.md5

    def calculate_size(self) -> None:
        """ Calculate the number of bytes in a given file.
        """
        self.size_in_bytes = calculate_file_size(self.absolute_file_path)
        return self.size_in_bytes

    def get_display_name(self):
        """ For dev convenience """
        if not self.filename:
            return self.source_filename
        else:
            return self.filename

    def get_extension(self):
        """ Returns the lowercased extension of the filename
        Thanks to https://stackoverflow.com/a/541408/763705 """
        return FileUtils.get_extension(self.filename)

    def is_blacklisted(self):
        return self.get_extension() in [".xml", ".chp", ".exp"]

    def delete_local_file(self):
        """ Deletes this file from the local file system."""
        try:
            os.remove(self.absolute_file_path)
        except OSError:
            pass
        except TypeError:
            pass
        except Exception:
            logger.exception(
                "Unexpected delete file exception.", absolute_file_path=self.absolute_file_path
            )
        self.is_downloaded = False
        self.save()

    def has_blocking_jobs(self, own_processor_id=None) -> bool:
        # Ignore transcriptome jobs because they always queue two jobs
        # for the same file.
        transcriptome_job_types = [
            ProcessorPipeline.TRANSCRIPTOME_INDEX_LONG.value,
            ProcessorPipeline.TRANSCRIPTOME_INDEX_SHORT.value,
        ]

        # If the file has a processor job that should not have been
        # retried, then it still shouldn't be retried.
        # Exclude the ones that were aborted.
        no_retry_processor_jobs = self.processor_jobs.filter(
            no_retry=True, worker_version=SYSTEM_VERSION
        ).exclude(abort=True, pipeline_applied__in=transcriptome_job_types)

        # If the file has a processor job that hasn't even started
        # yet, then it doesn't need another.
        incomplete_processor_jobs = self.processor_jobs.filter(
            end_time__isnull=True, success__isnull=True, retried=False
        ).exclude(pipeline_applied__in=transcriptome_job_types)

        if own_processor_id:
            incomplete_processor_jobs = incomplete_processor_jobs.exclude(id=own_processor_id)

        # Check if there's any jobs which should block another
        # processing attempt.
        blocking_jobs = no_retry_processor_jobs | incomplete_processor_jobs

        return blocking_jobs.first() is not None

    def needs_processing(self, own_processor_id=None) -> bool:
        """Returns False if original_file has been or is being processed.

        Returns True otherwise.

        If own_processor_id is supplied then it will be ignored so
        that processor jobs can use this function without their job
        being counted as currently processing this file.
        """
        if self.has_blocking_jobs(own_processor_id):
            return False

        sample = self.samples.first()
        if not sample:
            return True

        if sample.source_database == "SRA":
            computed_file = sample.get_most_recent_smashable_result_file()

            # If there's no smashable file then we should check the quant.sf file.
            if not computed_file:
                computed_file = sample.get_most_recent_quant_sf_file()

            # If there's neither a quant.sf file nor a smashable file
            # then we definitely need to process it.
            if not computed_file:
                return True

            if (
                computed_file.s3_bucket
                and computed_file.s3_key
                and computed_file.result.organism_index is not None
                and computed_file.result.organism_index.salmon_version == CURRENT_SALMON_VERSION
            ):
                # If the file wasn't computed with the latest
                # version of salmon, then it should be rerun
                # with the latest version of salmon.
                return False
        else:
            # If this original_file has multiple samples (is an
            # archive), and any of them haven't been processed, we'll
            # need the entire archive in order to process any of them.
            # A check to not re-processed the already processed
            # samples in the archive will happen elsewhere before
            # dispatching.
            for sample in self.samples.all():
                if not sample.is_processed:
                    return True
                computed_file = sample.get_most_recent_smashable_result_file()
                if not computed_file:
                    return True
                if settings.RUNNING_IN_CLOUD and (
                    computed_file.s3_bucket is None or computed_file.s3_key is None
                ):
                    return True

            return False

        # If we aren't sure, prefer reprocessing over never processing.
        return True

    def needs_downloading(self, own_processor_id=None) -> bool:
        """Determine if a file needs to be downloaded.

        This is true if the file has already been downloaded and lost
        without getting processed.
        """
        # If the file is downloaded and the file actually exists on disk,
        # then it doens't need to be downloaded.
        if self.absolute_file_path and os.path.exists(self.absolute_file_path):
            # ok a file exists, if this file has an SHA1 ensure that it's the same
            existing_file_sha1 = calculate_sha1(self.absolute_file_path)
            if self.sha1 and self.sha1 != existing_file_sha1:
                return True
            # otherwise, sha1 matches and the file doesn't need to be downloaded
            return False

        unstarted_downloader_jobs = self.downloader_jobs.filter(
            start_time__isnull=True, success__isnull=True, retried=False
        )

        # If the file has a downloader job that hasn't even started yet,
        # then it doesn't need another.
        if unstarted_downloader_jobs.count() > 0:
            return False

        # Do an extra check for blocking jobs for trancsriptome indices.
        # This is necessary because needs_processing() won't check
        # the blocking jobs for them because they're supposed to
        # have multiple processor jobs. However if the file does need to
        # be redownloaded, we only want one downloader job to be recreated.
        if self.has_blocking_jobs(own_processor_id):
            return False

        # If this file has been processed, then it doesn't need to be downloaded again.
        return self.needs_processing(own_processor_id)

    def is_affy_data(self) -> bool:
        """Return true if original_file is a CEL file or a gzipped CEL file.
        """
        upper_name = self.source_filename.upper()
        return (len(upper_name) > 4 and upper_name[-4:] == ".CEL") or (
            len(upper_name) > 7 and upper_name[-7:] == ".CEL.GZ"
        )

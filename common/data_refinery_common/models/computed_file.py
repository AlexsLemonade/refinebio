import os
import shutil

from django.conf import settings
from django.db import models
from django.utils import timezone

import boto3
from botocore.client import Config

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.managers import PublicObjectsManager
from data_refinery_common.utils import calculate_file_size, calculate_sha1

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))

logger = get_and_configure_logger(__name__)


class ComputedFile(models.Model):
    """ A representation of a file created by a data-refinery process """

    class Meta:
        db_table = "computed_files"
        get_latest_by = "created_at"

        indexes = [
            models.Index(fields=["filename"]),
        ]

    def __str__(self):
        return "ComputedFile: " + str(self.filename)

    SVD_ALGORITHM_CHOICES = (
        ("NONE", "None"),
        ("RANDOMIZED", "randomized"),
        ("ARPACK", "arpack"),
    )

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Object relations
    samples = models.ManyToManyField("Sample", through="SampleComputedFileAssociation")

    # File related
    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    # TODO: make this work w/ migrations:
    # absolute_file_path = models.CharField(max_length=255)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    # Relations
    result = models.ForeignKey(
        "ComputationalResult", blank=False, null=False, on_delete=models.CASCADE
    )

    # Scientific
    is_smashable = models.BooleanField(default=False)
    is_qc = models.BooleanField(default=False)
    is_qn_target = models.BooleanField(default=False)

    # Compendia details
    quant_sf_only = models.BooleanField(default=False)
    is_compendia = models.BooleanField(default=False)
    svd_algorithm = models.CharField(
        max_length=255,
        choices=SVD_ALGORITHM_CHOICES,
        default="NONE",
        help_text="The SVD algorithm that was used to generate the file.",
    )
    compendia_organism = models.ForeignKey(
        "Organism", blank=True, null=True, on_delete=models.CASCADE
    )
    compendia_version = models.IntegerField(blank=True, null=True)

    # AWS
    s3_bucket = models.CharField(max_length=255, blank=True, null=True)
    s3_key = models.CharField(max_length=255, blank=True, null=True)

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
        return super(ComputedFile, self).save(*args, **kwargs)

    def sync_to_s3(self, s3_bucket, s3_key) -> bool:
        """ Syncs a file to AWS S3.
        """
        if not settings.RUNNING_IN_CLOUD:
            return True

        try:
            S3.upload_file(
                self.absolute_file_path,
                s3_bucket,
                s3_key,
                ExtraArgs={"ACL": "public-read", "StorageClass": "STANDARD_IA"},
            )
        except Exception:
            logger.exception(
                "Error uploading computed file to S3",
                computed_file_id=self.pk,
                s3_key=self.s3_key,
                s3_bucket=self.s3_bucket,
            )
            return False

        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        self.save()

        return True

    def sync_from_s3(self, force=False, path=None):
        """ Downloads a file from S3 to the local file system.
        Returns the absolute file path.
        """
        path = path if path is not None else self.absolute_file_path

        if not settings.RUNNING_IN_CLOUD and not force:
            if os.path.exists(path):
                return path
            else:
                # If the file doesn't exist at path and we're not
                # running in the cloud, then the file is almost
                # certainly at its absolute_file_path because it never got deleted.
                if os.path.exists(self.absolute_file_path):
                    shutil.copyfile(self.absolute_file_path, path)
                    return path
                else:
                    # We don't have the file :(
                    return None

        target_directory = os.path.dirname(path)
        os.makedirs(target_directory, exist_ok=True)

        if not self.s3_bucket or not self.s3_key:
            raise ValueError("Tried to download a computed file with no s3_bucket or s3_key")

        try:
            S3.download_file(self.s3_bucket, self.s3_key, path)

            # Veryify sync integrity
            synced_sha1 = calculate_sha1(path)

            if self.sha1 != synced_sha1:
                raise AssertionError("SHA1 of downloaded ComputedFile doesn't match database SHA1!")

            return path
        except Exception as e:
            logger.exception(e, computed_file_id=self.pk)
            return None

    def change_s3_location(self, new_bucket: str, new_key: str) -> bool:
        """Moves the file from its current location in S3.

        The new location will be set based on `new_bucket` and
        `new_key`. The s3_bucket and s3_key properties will be updated
        to reflect this on a successful move.
        """
        old_bucket = self.s3_bucket
        old_key = self.s3_key
        copy_source = {"Bucket": old_bucket, "Key": old_key}
        try:
            response = S3.copy_object(Bucket=new_bucket, CopySource=copy_source, Key=new_key)
        except:
            logger.exception(
                "Could not copy computed file within S3",
                computed_file_id=self.id,
                source_bucket=old_bucket,
                source_key=old_key,
                destination_bucket=new_bucket,
                destination_key=new_key,
            )
            return False

        try:
            self.s3_bucket = new_bucket
            self.s3_key = new_key
            self.save()
        except:
            logger.exception(
                "Could not save computed file after it was copied!!!",
                computed_file_id=self.id,
                source_bucket=old_bucket,
                source_key=old_key,
                destination_bucket=new_bucket,
                destination_key=new_key,
            )
            return False

        try:
            response = S3.delete_object(Bucket=old_bucket, Key=old_key)
        except:
            logger.exception(
                "Could not delete computed file after it was copied and saved!!!",
                computed_file_id=self.id,
                source_bucket=old_bucket,
                source_key=old_key,
                destination_bucket=new_bucket,
                destination_key=new_key,
            )
            return False

        return True

    def calculate_sha1(self) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        self.sha1 = calculate_sha1(self.absolute_file_path)
        return self.sha1

    def calculate_size(self) -> None:
        """ Calculate the number of bytes in a given file.
        """
        self.size_in_bytes = calculate_file_size(self.absolute_file_path)
        return self.size_in_bytes

    def delete_local_file(self, force=False):
        """ Deletes a file from the path and actually removes it from the file system."""
        if not settings.RUNNING_IN_CLOUD and not force:
            return

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

    def delete_s3_file(self, force=False):
        # If we're not running in the cloud then we shouldn't try to
        # delete something from S3 unless force is set.
        if not settings.RUNNING_IN_CLOUD and not force:
            return False

        try:
            S3.delete_object(Bucket=self.s3_bucket, Key=self.s3_key)
        except:
            logger.exception(
                "Failed to delete S3 object for Computed File.",
                computed_file=self.id,
                s3_object=self.s3_key,
            )
            return False

        self.s3_key = None
        self.s3_bucket = None
        self.save()
        return True

    def get_synced_file_path(self, force=False, path=None):
        """ Fetches the absolute file path to this ComputedFile, fetching from S3 if it
        isn't already available locally. """
        if path:
            if os.path.exists(path):
                return path
            else:
                return self.sync_from_s3(force, path)
        else:
            if os.path.exists(self.absolute_file_path):
                return self.absolute_file_path
            else:
                return self.sync_from_s3(force)

    @property
    def s3_url(self):
        """ Render the resulting HTTPS URL for the S3 object."""
        return self.get_s3_url()

    def get_s3_url(self):
        """ Render the resulting HTTPS URL for the S3 object."""
        if (self.s3_key) and (self.s3_bucket):
            return "https://s3.amazonaws.com/" + self.s3_bucket + "/" + self.s3_key
        else:
            return None

    @property
    def download_url(self):
        """ A temporary URL from which the file can be downloaded. """
        return self.create_download_url()

    def create_download_url(self):
        """ Create a temporary URL from which the file can be downloaded."""
        if settings.RUNNING_IN_CLOUD and self.s3_bucket and self.s3_key:
            return S3.generate_presigned_url(
                ClientMethod="get_object",
                Params={"Bucket": self.s3_bucket, "Key": self.s3_key},
                ExpiresIn=(60 * 60 * 7 * 24),  # 7 days in seconds.
            )
        else:
            return None

    def has_been_log2scaled(self):
        """ Return true if this is a smashable file that has been log2 scaled """
        return self.is_smashable and self.filename.endswith("lengthScaledTPM.tsv")

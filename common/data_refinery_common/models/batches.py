from enum import Enum
import os
import urllib
import shutil
import boto3
from django.db import transaction
from django.db import models
from data_refinery_common.models.base_models import TimeTrackedModel
from data_refinery_common.models.surveys import SurveyJob
from data_refinery_common.utils import get_env_variable

# Import and set logger
import logging
logger = logging.getLogger(__name__)


DEFAULT_BATCH_PREFIX = "batch_"
RAW_PREFIX = get_env_variable("RAW_PREFIX")
TEMP_PREFIX = get_env_variable("TEMP_PREFIX")
PROCESSED_PREFIX = get_env_variable("PROCESSED_PREFIX")
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR")
USE_S3 = get_env_variable("USE_S3") == "True"
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME")


class BatchStatuses(Enum):
    """Valid values for the status field of the Batch model."""

    NEW = "NEW"
    DOWNLOADED = "DOWNLOADED"
    PROCESSED = "PROCESSED"


class Batch(TimeTrackedModel):
    """Represents a batch of data.

    The definition of a Batch is intentionally that vague. What a batch
    is will vary from source to source. It could be a single file, or
    a group of files with some kind of logical grouping such as an
    experiment.
    """

    survey_job = models.ForeignKey(SurveyJob, on_delete=models.PROTECT)
    source_type = models.CharField(max_length=256)
    pipeline_required = models.CharField(max_length=256)
    platform_accession_code = models.CharField(max_length=32)
    # One source type uses organism names for this field and the
    # longest organism name is apparently Parastratiosphecomyia
    # stratiosphecomyioides at 44 characters.
    experiment_accession_code = models.CharField(max_length=64)
    experiment_title = models.CharField(max_length=256)
    status = models.CharField(max_length=20)
    release_date = models.DateField()
    last_uploaded_date = models.DateField()

    # This corresponds to the organism taxonomy ID from NCBI.
    organism_id = models.IntegerField()
    # This is the organism name as it appeared in the experiment.
    organism_name = models.CharField(max_length=256)

    _files = None

    def _get_files(self):
        if self._files is None:
            return list(self.file_set.all())
        else:
            return self._files

    def _set_files(self, new_files):
        self._files = new_files

    files = property(_get_files, _set_files)

    @classmethod
    def is_new_batch(cls, batch) -> bool:
        file_names = batch.file_set.all().values("name")
        matching_batches = cls.objects.filter(name__in=file_names)
        return matching_batches.count() == 0

    @classmethod
    @transaction.atomic
    def create_file_and_relationships(cls, *args, **kwargs):
        """Inits and saves a batch and its relationships to its files.

        Expects keyword arguments that could be passed to the init
        method for Batch, with the addition of the keyword
        argument 'files', which must be specified as a list of File
        objects which have already been saved (so they have an id).
        """
        files = kwargs.pop('files', None)
        this_batch = cls(*args, **kwargs)
        if files is None:
            raise KeyError("The 'files' argument must be specified.")
        else:
            this_batch.save()
            this_batch.files.add(*files)

        return this_batch

    class Meta:
        db_table = "batches"


class File(TimeTrackedModel):
    name = models.CharField(max_length=1024)
    download_url = models.CharField(max_length=4096)
    raw_format = models.CharField(max_length=256, null=True)
    processed_format = models.CharField(max_length=256, null=True)
    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    size_in_bytes = models.BigIntegerField()

    # This field will denote where in our system the file can be found.
    internal_location = models.CharField(max_length=256, null=True)

    def _get_downloaded_name(self) -> str:
        path = urllib.parse.urlparse(self.download_url).path
        return os.path.basename(path)

    def get_base_name(self) -> str:
        return self.name.replace(("." + self.raw_format), "")

    def get_processed_name(self) -> str:
        file_base = self.get_base_name()
        return file_base + "." + self.processed_format

    def get_raw_dir(self) -> str:
        if USE_S3:
            return os.path.join(RAW_PREFIX, self.internal_location)
        else:
            return os.path.join(LOCAL_ROOT_DIR, RAW_PREFIX, self.internal_location)

    def get_raw_download_path(self) -> str:
        """Get the path to the downloaded file in the raw directory.

        In cases where extraction is necessary, this will not match the
        name of the batch's extracted file.
        """
        return os.path.join(self.get_raw_dir(), self._get_downloaded_name())

    def get_raw_path(self) -> str:
        return os.path.join(self.get_raw_dir(), self.name)

    # The following four functions use the ID of the batch in the
    # temporary paths so it can be removed after processing is complete
    # without interfering with other jobs. An optional dir_name can be
    # used instead.
    def get_temp_dir(self, dir_name: str=None) -> str:
        dir_name = dir_name if dir_name is not None else DEFAULT_BATCH_PREFIX + str(self.batch.id)
        return os.path.join(LOCAL_ROOT_DIR,
                            TEMP_PREFIX,
                            self.internal_location,
                            dir_name)

    def get_temp_download_path(self, dir_name: str=None) -> str:
        """Returns the path of the downloaded file in the temp directory.

        The base name of returned path is set to the base name of
        download_url. In cases where extraction is necessary, this
        will not match the name of the batch's extracted file.
        """
        return os.path.join(self.get_temp_dir(dir_name),
                            self._get_downloaded_name())

    def get_temp_pre_path(self, dir_name: str=None) -> str:
        """Returns the path of the pre-processed file for the batch."""
        return os.path.join(self.get_temp_dir(dir_name),
                            self.name)

    def get_temp_post_path(self, dir_name: str=None) -> str:
        """Returns the path of the post-processed file for the batch."""
        return os.path.join(self.get_temp_dir(dir_name),
                            self.get_processed_name())

    def get_processed_dir(self) -> str:
        if USE_S3:
            return os.path.join(PROCESSED_PREFIX, self.internal_location)
        else:
            return os.path.join(LOCAL_ROOT_DIR, PROCESSED_PREFIX, self.internal_location)

    def get_processed_path(self) -> str:
        return os.path.join(self.get_processed_dir(),
                            self.get_processed_name())

    def _upload_file(self, from_path: str, to_dir: str, to_path: str) -> None:
        """Move the file from from_path to to_path.

        Depending on the value of the USE_S3 environment variable this
        will either move the file to a local directory or to S3.
        """
        if USE_S3:
            bucket = boto3.resource("s3").Bucket(S3_BUCKET_NAME)
            with open(from_path, 'rb') as from_file:
                bucket.put_object(Key=to_path, Body=from_file)
        else:
            os.makedirs(to_dir, exist_ok=True)
            shutil.copyfile(from_path, to_path)

    def upload_raw_file(self, dir_name: str=None) -> None:
        """Moves the batch's raw file out of the temp directory.

        Depending on the value of the USE_S3 environment variable this
        will either move the batch's raw file to the RAW_PREFIX directory
        or to S3.
        """
        temp_path = self.get_temp_pre_path(dir_name)
        raw_dir = self.get_raw_dir()
        raw_path = self.get_raw_path()

        logger.debug("Moving file from %s to %s.", temp_path, raw_path)
        self._upload_file(temp_path, raw_dir, raw_path)

    def download_raw_file(self, dir_name: str=None) -> str:
        """Moves the batch's raw file to the temp directory.

        Depending on the value of the USE_S3 environment variable this
        will either move the batch's raw file from the RAW_PREFIX directory
        or from S3.
        Returns the path the file was downloaded to.
        """
        raw_path = self.get_raw_path()
        temp_dir = self.get_temp_dir(dir_name)
        os.makedirs(temp_dir, exist_ok=True)
        temp_path = self.get_temp_pre_path(dir_name)
        if USE_S3:
            bucket = boto3.resource("s3").Bucket(S3_BUCKET_NAME)
            with open(temp_path, 'wb') as temp_file:
                bucket.download_fileobj(raw_path, temp_file)
        else:
            shutil.copyfile(raw_path, temp_path)

        return temp_path

    def download_processed_file(self, dir_name: str=None) -> str:
        """Moves the batch's processed file to the temp directory.

        Depending on the value of the USE_S3 environment variable this
        will either move the batch's processed file from the
        PROCESSED_PREFIX directory or from S3.
        Returns the path the file was downloaded to.
        """
        processed_path = self.get_processed_path()
        temp_dir = self.get_temp_dir(dir_name)
        os.makedirs(temp_dir, exist_ok=True)
        temp_path = self.get_temp_post_path(dir_name)
        if USE_S3:
            bucket = boto3.resource("s3").Bucket(S3_BUCKET_NAME)
            with open(temp_path, 'wb') as temp_file:
                bucket.download_fileobj(processed_path, temp_file)
        else:
            shutil.copyfile(processed_path, temp_path)

        return temp_path

    def upload_processed_file(self, dir_name: str=None) -> None:
        """Moves the batch's processed file out of the temp directory.

        Depending on the value of the USE_S3 environment variable this may
        just be to the PROCESSED_PREFIX directory or it may be to S3.
        """
        temp_path = self.get_temp_post_path(dir_name)
        processed_dir = self.get_processed_dir()
        processed_path = self.get_processed_path()
        self._upload_file(temp_path, processed_dir, processed_path)

    def remove_temp_directory(self, dir_name: str=None) -> None:
        temp_dir = self.get_temp_dir(dir_name)
        if os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir)

    def remove_raw_files(self) -> None:
        """Cleans up the raw files that were downloaded for a batch.

        Depending on the value of the USE_S3 environment variable the
        files cleaned up may be stored locally or within S3.
        """
        raw_path = self.get_raw_path()
        if USE_S3:
            bucket = boto3.resource("s3").Bucket(S3_BUCKET_NAME)
            bucket.delete_objects(
                Delete={
                    'Objects': [
                        {
                            'Key': raw_path
                        }
                    ]
                }
            )
        else:
            os.remove(raw_path)

    class Meta:
        db_table = "files"


class BatchKeyValue(TimeTrackedModel):
    """Tracks additional fields for Batches.

    Useful for fields that would be sparsely populated if they were
    their own columns. I.e. one source may have an extra field or two
    that are worth tracking but are specific to that source.
    """

    batch = models.ForeignKey(Batch, on_delete=models.CASCADE)
    key = models.CharField(max_length=256)
    value = models.CharField(max_length=256)

    class Meta:
        db_table = "batch_key_values"

import hashlib
import io
import os
import pytz

from datetime import datetime
from functools import partial

from django.contrib.postgres.fields import HStoreField
from django.db import transaction
from django.db import models
from django.utils import timezone

from data_refinery_common.models.organism import Organism

"""
# First Order Classes

This represent the primary data types we will be querying
and filtering against.

"""

class Sample(models.Model):
    """
    An individual sample.
    """

    class Meta:
        db_table = "samples"

    def __str__ (self):
        return "Sample: " + self.accession_code

    # Identifiers
    accession_code = models.CharField(max_length=255, unique=True)

    # Relations
    organism = models.ForeignKey(Organism, blank=True, null=True, on_delete=models.SET_NULL)

    # Historical Properties
    source_archive_url = models.CharField(max_length=255)
    source_filename = models.CharField(max_length=255, blank=False)
    source_absolute_file_path = models.CharField(max_length=255)

    # Scientific Properties
    has_raw = models.BooleanField(default=True) # Did this sample have a raw data source?

    # Crunch Properties
    is_downloaded = models.BooleanField(default=False)
    is_processed = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)


class SampleAnnotation(models.Model):
    """ Semi-standard information associated with a Sample """

    class Meta:
        db_table = "sample_annotations"

    # Relations
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class Experiment(models.Model):
    """ An Experiment or Study """

    class Meta:
        db_table = "experiments"

    def __str__ (self):
        return "Experiment: " + self.accession_code

    # Identifiers
    accession_code = models.CharField(max_length=64, unique=True)

    # Historical Properties
    source_database = models.CharField(max_length=32) # "ArrayExpress, "SRA"
    source_url = models.CharField(max_length=256)

    # Properties
    ## I was always under the impression that TextFields were slower
    ## than CharFields, however the Postgres documentation disagrees: 
    ## https://www.postgresql.org/docs/9.0/static/datatype-character.html
    title = models.TextField()
    description = models.TextField()
    protocol_description = models.TextField()
    platform_accession_code = models.CharField(max_length=256, blank=True)
    platform_name = models.CharField(max_length=256, blank=True)
    submitter_institution = models.CharField(max_length=256, blank=True)
    has_publication = models.BooleanField(default=False)
    publication_title = models.TextField()
    publication_doi = models.CharField(max_length=64, blank=True)
    pubmed_id = models.CharField(max_length=16, blank=True)
    source_first_published = models.CharField(max_length=256, blank=True)
    source_last_updated = models.CharField(max_length=256, blank=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class ExperimentAnnotation(models.Model):
    """ Semi-standard information associated with an Experiment """

    class Meta:
        db_table = "experiment_annotations"

    # Relations
    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class ComputationalResult(models.Model):
    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"

    def __str__ (self):
        return "ComputationalResult: " + str(self.pk)

    command_executed = models.CharField(max_length=255, blank=True)
    program_version = models.CharField(max_length=255)
    system_version = models.CharField(max_length=255) # Generally defined in from data_refinery_workers._version import __version__
    is_ccdl = models.BooleanField(default=True)

    # Stats
    time_start = models.DateTimeField(blank=True, null=True) 
    time_end = models.DateTimeField(blank=True, null=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now) 

class CompultationalResultAnnotation(models.Model):
    """ Non-standard information associated with an ComputationalResult """

    class Meta:
        db_table = "computational_result_annotations"

    # Relations
    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now) 

# TODO
# class Gene(models.Model):
    """ A representation of a Gene """

#     class Meta:
#         db_table = "genes"

class Index(models.Model):
    """ A special type of process result, necessary for processing other SRA samples """

    class Meta:
        db_table = "index"

    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)
    index_type = models.CharField(max_length=255) # XXX ex "TRANSCRIPTOME_LONG", "TRANSCRIPTOME_SHORT", ???
    source_version = models.CharField(max_length=255) # Where do we get this from
    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Common Properties
    is_public = models.BooleanField(default=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now) 

"""
# Files

These are the database representations of files
which live on local disk, on ephemeral storage,
or on AWS cloud services.
"""

class OriginalFile(models.Model):
    """ A representation of a file from an external source """

    class Meta:
        db_table = "original_files"

    def __str__ (self):
        return "OriginalFile: " + self.get_display_name()

    # XXX: move to `filename`
    file_name = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField(blank=True, null=True)
    sha1 = models.CharField(max_length=64)

    # Reference properties
    sample = models.ForeignKey(Sample, blank=True, null=True, on_delete=models.CASCADE)

    # Historical Properties
    source_url = models.CharField(max_length=255)
    is_archive = models.BooleanField(default=True)
    source_filename = models.CharField(max_length=255, blank=False)

    # Scientific Properties
    has_raw = models.BooleanField(default=True) # Did this sample have a raw data source?

    # Crunch Properties
    is_downloaded = models.BooleanField(default=False)
    is_processed = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now) 

    def calculate_sha1(self, absolute_file_path:str=None) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        if not absolute_file_path:
            absolute_file_path = self.absolute_file_path

        hash_object = hashlib.sha1() 
        with open(absolute_file_path, mode='rb') as open_file:
            for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b''):
                hash_object.update(buf)
        
        self.sha1 = hash_object.hexdigest()
        return self.sha1

    def calculate_size(self, absolute_file_path:str=None) -> None:
        """ Calculate the number of bites in a given file.
        """
        if not absolute_file_path:
            absolute_file_path = self.absolute_file_path

        self.size_in_bytes = os.path.getsize(absolute_file_path)
        return self.size_in_bytes

    def get_display_name(self):
        """ For dev convenience """
        if not self.file_name:
            return self.source_filename
        else:
            return self.file_name

class ComputedFile(models.Model):
    """ A representation of a file created by a data-refinery process """

    class Meta:
        db_table = "computed_files"

    def __str__ (self):
        return "ComputedFile: " + str(self.file_name)

    file_name = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)
    s3_bucket = models.CharField(max_length=255)
    s3_key = models.CharField(max_length=255)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now) 

    def sync_to_s3(self, s3_bucket=None, s3_key=None) -> bool:
        """ Syncs a file to AWS S3. 

        XXX: TODO!
        """
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        return True

    def calculate_sha1(self, absolute_file_path:str=None) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        if not absolute_file_path:
            absolute_file_path = self.absolute_file_path

        hash_object = hashlib.sha1() 
        with open(absolute_file_path, mode='rb') as open_file:
            for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b''):
                hash_object.update(buf)
        
        self.sha1 = hash_object.hexdigest()
        return self.sha1

    def calculate_size(self, absolute_file_path:str=None) -> None:
        """ Calculate the number of bites in a given file.
        """
        if not absolute_file_path:
            absolute_file_path = self.absolute_file_path

        self.size_in_bytes = os.path.getsize(absolute_file_path)
        return self.size_in_bytes

"""
# Associations

These represent the relationships between items in the other tables. 
"""

class ExperimentSampleAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_sample_associations"

class DownloaderJobOriginalFileAssociation(models.Model):

    downloader_job = models.ForeignKey("data_refinery_common.DownloaderJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "downloaderjob_originalfile_associations"

class ProcessorJobOriginalFileAssociation(models.Model):

    processor_job = models.ForeignKey("data_refinery_common.ProcessorJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "processorjob_originalfile_associations"

class SampleResultAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "sample_result_associations"

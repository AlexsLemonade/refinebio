import hashlib
import io
import os
import pytz

from datetime import datetime
from django.db import transaction
from django.db import models
from django.utils import timezone
from functools import partial

from data_refinery_common.models.organism import Organism

"""
# Miscellaneous

Everything that's left over.

"""

# class Organism(models.Model):
#     """Provides a lookup between organism name and taxonomy ids.

#     Should only be used via the two class methods get_name_for_id and
#     get_id_for_name. These methods will populate the database table
#     with any missing values by accessing the NCBI API.
#     """

#     taxonomy_id = models.IntegerField()
#     name = models.CharField(max_length=256)
#     is_scientific_name = models.BooleanField(default=False)

#     # @classmethod
#     # def get_name_for_id(cls, taxonomy_id: int) -> str:
#     #     try:
#     #         organism = (cls.objects
#     #                     .filter(taxonomy_id=taxonomy_id)
#     #                     .order_by("-is_scientific_name")
#     #                     [0])
#     #     except IndexError:
#     #         name = get_scientific_name(taxonomy_id).upper()
#     #         organism = Organism(name=name,
#     #                             taxonomy_id=taxonomy_id,
#     #                             is_scientific_name=True)
#     #         organism.save()
#     #
#     #     return organism.name
#     #
#     # @classmethod
#     # def get_id_for_name(cls, name: str) -> id:
#     #     name = name.upper()
#     #     try:
#     #         organism = (cls.objects
#     #                     .filter(name=name)
#     #                     [0])
#     #     except IndexError:
#     #         is_scientific_name = False
#     #         try:
#     #             taxonomy_id = get_taxonomy_id_scientific(name)
#     #             is_scientific_name = True
#     #         except UnscientificNameError:
#     #             taxonomy_id = get_taxonomy_id(name)
#     #
#     #         organism = Organism(name=name,
#     #                             taxonomy_id=taxonomy_id,
#     #                             is_scientific_name=is_scientific_name)
#     #         organism.save()
#     #
#     #     return organism.taxonomy_id

#     class Meta:
#         db_table = "organisms"

"""
# First Order Classes

This represent the primary data types we will be querying
and filtering against.

"""

class Sample(models.Model):
    """
    An individual sample.
    Examples:

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

    class Meta:
        db_table = "sample_annotations"

    # Relations
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    key = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class Experiment(models.Model):

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
    title = models.CharField(max_length=256)
    description = models.TextField()
    platform_accession_code = models.CharField(max_length=256)
    platform_name = models.CharField(max_length=256)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class ExperimentAnnotation(models.Model):

    class Meta:
        db_table = "experiment_annotations"

    # Relations
    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    key = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

class ComputationalResult(models.Model):

    class Meta:
        db_table = "computational_results"

    def __str__ (self):
        return "ComputationalResult: " + str(self.pk)

    command_executed = models.CharField(max_length=255, blank=True)
    program_version = models.CharField(max_length=255) # Define in settings!
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

    class Meta:
        db_table = "computational_result_annotations"

    # Relations
    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    key = models.CharField(max_length=255)
    value = models.CharField(max_length=255)

class Gene(models.Model):

    class Meta:
        db_table = "genes"

class Index(models.Model):

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

    class Meta:
        db_table = "computed_files"

    def __str__ (self):
        return "ComputedFile: " + str(self.file_name)

    file_name = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)
    s3_bucket = models.CharField(max_length=64)
    s3_key = models.CharField(max_length=64)

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


# XXX: Can't we just make lots of SampleResultAssociations instead?
# class ExperimentResultAssociation(models.Model):

#     class Meta:
#         db_table = "experiment_result_associations"


"""
# Jobs

These models correspond to the status of various task workers
responsible for the downloading and transformation of our
input data.

"""



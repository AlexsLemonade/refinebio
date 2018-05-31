import hashlib
import io
import os
import pytz
import uuid

from datetime import datetime
from functools import partial

from django.contrib.postgres.fields import HStoreField, JSONField
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

    def __str__(self):
        return "Sample: " + self.accession_code

    # Identifiers
    accession_code = models.CharField(max_length=255, unique=True)
    title = models.CharField(max_length=255, unique=False, blank=True)

    # Relations
    organism = models.ForeignKey(Organism, blank=True, null=True, on_delete=models.SET_NULL)
    results = models.ManyToManyField('ComputationalResult', through='SampleResultAssociation')
    original_files = models.ManyToManyField('OriginalFile', through='OriginalFileSampleAssociation')

    # Historical Properties
    source_database = models.CharField(max_length=255, blank=False)
    source_archive_url = models.CharField(max_length=255)
    source_filename = models.CharField(max_length=255, blank=False)
    source_absolute_file_path = models.CharField(max_length=255)
    has_raw = models.BooleanField(default=True)  # Did this sample have a raw data source?

    # Technological Properties
    platform_accession_code = models.CharField(max_length=256, blank=True)
    platform_name = models.CharField(max_length=256, blank=True)
    technology = models.CharField(max_length=256, blank=True)

    # Scientific Properties
    sex = models.CharField(max_length=255, blank=True)
    age = models.DecimalField(max_length=255, blank=True, max_digits=8, decimal_places=3, null=True)
    specimen_part = models.CharField(max_length=255, blank=True)
    genotype = models.CharField(max_length=255, blank=True)
    disease = models.CharField(max_length=255, blank=True)
    disease_stage = models.CharField(max_length=255, blank=True)
    cell_line = models.CharField(max_length=255, blank=True)
    treatment = models.CharField(max_length=255, blank=True)
    race = models.CharField(max_length=255, blank=True)
    subject = models.CharField(max_length=255, blank=True)
    compound = models.CharField(max_length=255, blank=True)
    time = models.CharField(max_length=255, blank=True)

    # Crunch Properties
    is_downloaded = models.BooleanField(default=False)
    is_processed = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(Sample, self).save(*args, **kwargs)


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

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(SampleAnnotation, self).save(*args, **kwargs)


class Experiment(models.Model):

    """ An Experiment or Study """

    class Meta:
        db_table = "experiments"

    def __str__(self):
        return "Experiment: " + self.accession_code

    # Relations
    samples = models.ManyToManyField('Sample', through='ExperimentSampleAssociation')
    organisms = models.ManyToManyField('Organism', through='ExperimentOrganismAssociation')

    # Identifiers
    accession_code = models.CharField(max_length=64, unique=True)

    # Historical Properties
    source_database = models.CharField(max_length=32)  # "ArrayExpress, "SRA"
    source_url = models.CharField(max_length=256)

    # Properties
    ## I was always under the impression that TextFields were slower
    ## than CharFields, however the Postgres documentation disagrees:
    ## https://www.postgresql.org/docs/9.0/static/datatype-character.html
    title = models.TextField()
    description = models.TextField()
    protocol_description = models.TextField(default="")
    technology = models.CharField(max_length=256, blank=True)
    submitter_institution = models.CharField(max_length=256, blank=True)
    has_publication = models.BooleanField(default=False)
    publication_title = models.TextField(default="")
    publication_doi = models.CharField(max_length=64, blank=True)
    pubmed_id = models.CharField(max_length=32, blank=True)
    source_first_published = models.DateTimeField(null=True)
    source_last_modified = models.DateTimeField(null=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(Experiment, self).save(*args, **kwargs)


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

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ExperimentAnnotation, self).save(*args, **kwargs)


class ComputationalResult(models.Model):

    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"

    def __str__(self):
        return "ComputationalResult: " + str(self.pk)

    command_executed = models.TextField(blank=True)
    program_version = models.TextField(blank=True)
    system_version = models.CharField(
        max_length=255)  # Generally defined in from data_refinery_workers._version import __version__
    is_ccdl = models.BooleanField(default=True)

    # Stats
    time_start = models.DateTimeField(blank=True, null=True)
    time_end = models.DateTimeField(blank=True, null=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ComputationalResult, self).save(*args, **kwargs)


class ComputationalResultAnnotation(models.Model):

    """ Non-standard information associated with an ComputationalResult """

    class Meta:
        db_table = "computational_result_annotations"

    # Relations
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=True)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ComputationalResultAnnotation, self).save(*args, **kwargs)

# TODO
# class Gene(models.Model):
    """ A representation of a Gene """

#     class Meta:
#         db_table = "genes"

class OrganismIndex(models.Model):
    """ A special type of process result, necessary for processing other SRA samples """

    class Meta:
        db_table = "organism_index"

    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)
    index_type = models.CharField(max_length=255) # ex., "TRANSCRIPTOME_LONG", "TRANSCRIPTOME_SHORT"
    source_version = models.CharField(max_length=255) # Where do we get this from
    result = models.ForeignKey(ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

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
        return super(OrganismIndex, self).save(*args, **kwargs)

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

    def __str__(self):
        return "OriginalFile: " + self.get_display_name()

    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField(blank=True, null=True)
    sha1 = models.CharField(max_length=64)

    # Relations
    samples = models.ManyToManyField('Sample', through='OriginalFileSampleAssociation')

    # Historical Properties
    source_url = models.CharField(max_length=255)
    is_archive = models.BooleanField(default=True)
    source_filename = models.CharField(max_length=255, blank=False)

    # Scientific Properties
    has_raw = models.BooleanField(default=True)  # Did this sample have a raw data source?

    # Crunch Properties
    is_downloaded = models.BooleanField(default=False)
    is_processed = models.BooleanField(default=False)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(OriginalFile, self).save(*args, **kwargs)

    def calculate_sha1(self) -> None:
        """ Calculate the SHA1 value of a given file.
        """

        hash_object = hashlib.sha1()
        with open(self.absolute_file_path, mode='rb') as open_file:
            for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b''):
                hash_object.update(buf)

        self.sha1 = hash_object.hexdigest()
        return self.sha1

    def calculate_size(self) -> None:
        """ Calculate the number of bytes in a given file.
        """
        self.size_in_bytes = os.path.getsize(self.absolute_file_path)
        return self.size_in_bytes

    def get_display_name(self):
        """ For dev convenience """
        if not self.filename:
            return self.source_filename
        else:
            return self.filename


class ComputedFile(models.Model):

    """ A representation of a file created by a data-refinery process """

    class Meta:
        db_table = "computed_files"

    def __str__(self):
        return "ComputedFile: " + str(self.filename)

    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)
    s3_bucket = models.CharField(max_length=255)
    s3_key = models.CharField(max_length=255)

    # Common Properties
    is_public = models.BooleanField(default=False)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(ComputedFile, self).save(*args, **kwargs)

    def sync_to_s3(self, s3_bucket=None, s3_key=None) -> bool:
        """ Syncs a file to AWS S3.

        XXX: TODO!
        """
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        return True

    def calculate_sha1(self) -> None:
        """ Calculate the SHA1 value of a given file.
        """
        hash_object = hashlib.sha1()
        with open(self.absolute_file_path, mode='rb') as open_file:
            for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b''):
                hash_object.update(buf)

        self.sha1 = hash_object.hexdigest()
        return self.sha1

    def calculate_size(self) -> None:
        """ Calculate the number of bytes in a given file.
        """
        self.size_in_bytes = os.path.getsize(self.absolute_file_path)
        return self.size_in_bytes


class Dataset(models.Model):

    """ A Dataset is a desired set of experiments/samples to smash and download """

    AGGREGATE_CHOICES = (
        ('EXPERIMENT', 'Experiment'),
        ('SPECIES', 'Species')
    )

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Experiments and samples live here: {'E-ABC-1': ['SAMP1', 'SAMP2']}
    # This isn't going to be queryable, so we can use JSON-in-text, just make
    # sure we validate properly in and out!
    data = JSONField(default={})

    # Processing properties
    aggregate_by = models.CharField(max_length=255, choices=AGGREGATE_CHOICES, default="EXPERIMENT")

    # State properties
    is_processing = models.BooleanField(default=False)  # Data is still editable
    is_processed = models.BooleanField(default=False)  # Result has been made
    is_available = models.BooleanField(default=False)  # Result is ready for delivery

    # Delivery properties
    email_address = models.CharField(max_length=255, blank=True, null=True)
    expires_on = models.DateTimeField(blank=True, null=True)

    # Deliverables
    s3_bucket = models.CharField(max_length=255)
    s3_key = models.CharField(max_length=255)

    # Common Properties
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(Dataset, self).save(*args, **kwargs)

"""
# Associations

These represent the relationships between items in the other tables.
"""


class ExperimentSampleAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_sample_associations"


class ExperimentOrganismAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_organism_associations"


class DownloaderJobOriginalFileAssociation(models.Model):

    downloader_job = models.ForeignKey(
        "data_refinery_common.DownloaderJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "downloaderjob_originalfile_associations"


class ProcessorJobOriginalFileAssociation(models.Model):

    processor_job = models.ForeignKey(
        "data_refinery_common.ProcessorJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "processorjob_originalfile_associations"


class OriginalFileSampleAssociation(models.Model):

    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "original_file_sample_associations"


class SampleResultAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "sample_result_associations"

import hashlib
import io
import os
import pytz
import uuid
import boto3

from botocore.client import Config

from datetime import datetime
from functools import partial

from django.conf import settings
from django.contrib.postgres.fields import ArrayField, HStoreField, JSONField
from django.db import transaction
from django.db import models
from django.utils import timezone

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.organism import Organism
from data_refinery_common.utils import get_env_variable, get_s3_url

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client('s3', config=Config(signature_version='s3v4'))

logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
CHUNK_SIZE = 1024 * 256 # chunk_size is in bytes

"""
# First Order Classes

This represent the primary data types we will be querying
and filtering against.

"""

class PublicObjectsManager(models.Manager):
    """
    Only returns objects that have is_public
    """
    def get_queryset(self):
        return super().get_queryset().filter(is_public=True)

class ProcessedObjectsManager(models.Manager):
    """
    Only returns objects that have is_processed and is_public
    """
    def get_queryset(self):
        return super().get_queryset().filter(is_processed=True, is_public=True)

class Sample(models.Model):
    """
    An individual sample.
    """

    class Meta:
        db_table = "samples"
        base_manager_name = 'public_objects'

    def __str__(self):
        return self.accession_code

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()
    processed_objects = ProcessedObjectsManager()

    # Identifiers
    accession_code = models.CharField(max_length=255, unique=True)
    title = models.CharField(max_length=255, unique=False, blank=True)

    # Relations
    organism = models.ForeignKey(Organism, blank=True, null=True, on_delete=models.SET_NULL)
    results = models.ManyToManyField('ComputationalResult', through='SampleResultAssociation')
    original_files = models.ManyToManyField('OriginalFile', through='OriginalFileSampleAssociation')
    computed_files = models.ManyToManyField('ComputedFile', through='SampleComputedFileAssociation')

    # Historical Properties
    source_database = models.CharField(max_length=255, blank=False)
    source_archive_url = models.CharField(max_length=255)
    source_filename = models.CharField(max_length=255, blank=False)
    source_absolute_file_path = models.CharField(max_length=255)
    has_raw = models.BooleanField(default=True)  # Did this sample have a raw data source?

    # Technological Properties
    platform_accession_code = models.CharField(max_length=256, blank=True)
    platform_name = models.CharField(max_length=256, blank=True)
    technology = models.CharField(max_length=256, blank=True) # MICROARRAY, RNA-SEQ
    manufacturer = models.CharField(max_length=256, blank=True)

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
    # TODO: rm is_downloaded from sample, it is on original_file instead
    is_downloaded = models.BooleanField(default=False)
    is_processed = models.BooleanField(default=False)

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
        return super(Sample, self).save(*args, **kwargs)

    def to_metadata_dict(self):
        """ Render this Sample as a dict """
        metadata = {}
        metadata['title'] = self.title
        metadata['accession_code'] = self.accession_code
        metadata['organism'] = self.organism.name
        metadata['source_archive_url'] = self.source_archive_url
        metadata['sex'] = self.sex
        metadata['age'] = self.age or ''
        metadata['specimen_part'] = self.specimen_part
        metadata['genotype'] = self.genotype
        metadata['disease'] = self.disease
        metadata['disease_stage'] = self.disease_stage
        metadata['cell_line'] = self.cell_line
        metadata['treatment'] = self.treatment
        metadata['race'] = self.race
        metadata['subject'] = self.subject
        metadata['compound'] = self.compound
        metadata['time'] = self.time
        metadata['platform'] = self.pretty_platform
        metadata['annotations'] = [data for data in self.sampleannotation_set.all().values_list('data', flat=True)]

        return metadata

    def get_result_files(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return self.computed_files.all()

    def get_most_recent_smashable_result_file(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return self.computed_files.filter(
                        is_smashable=True
                    ).first()

    @property
    def pipelines(self):
        """ Returns a list of related pipelines """
        return [p for p in self.results.values_list('pipeline', flat=True).distinct()]

    @property
    def pretty_platform(self):
        """ Turns

        [HT_HG-U133_Plus_PM] Affymetrix HT HG-U133+ PM Array Plate

        into

        Affymetrix HT HG-U133+ PM Array Plate (hthgu133pluspm)

        """
        if ']' in self.platform_name:
            platform_base = self.platform_name.split(']')[1].strip()
        else:
            platform_base = self.platform_name
        return platform_base + ' (' + self.platform_accession_code + ')'

class SampleAnnotation(models.Model):
    """ Semi-standard information associated with a Sample """

    class Meta:
        db_table = "sample_annotations"
        base_manager_name = 'public_objects'

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=False)

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
        return super(SampleAnnotation, self).save(*args, **kwargs)


class ProcessedPublicObjectsManager(models.Manager):
    """
    Only returns Experiments that are is_public and have related is_processed Samples.
    """
    def get_queryset(self):
        return super().get_queryset().filter(
            is_public=True,
            samples__is_processed=True,
            samples__is_public=True).distinct()


class Experiment(models.Model):
    """ An Experiment or Study """

    class Meta:
        db_table = "experiments"
        base_manager_name = 'public_objects'

    def __str__(self):
        return "Experiment: " + self.accession_code

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()
    processed_public_objects = ProcessedPublicObjectsManager()

    # Relations
    samples = models.ManyToManyField('Sample', through='ExperimentSampleAssociation')
    organisms = models.ManyToManyField('Organism', through='ExperimentOrganismAssociation')

    # Identifiers
    accession_code = models.CharField(max_length=64, unique=True)

    # Historical Properties
    source_database = models.CharField(max_length=32)  # "ArrayExpress, "SRA"
    source_url = models.CharField(max_length=256)

    # Properties
    # I was always under the impression that TextFields were slower
    # than CharFields, however the Postgres documentation disagrees:
    # https://www.postgresql.org/docs/9.0/static/datatype-character.html
    title = models.TextField()
    description = models.TextField()
    protocol_description = models.TextField(default="")
    technology = models.CharField(max_length=256, blank=True)
    submitter_institution = models.CharField(max_length=256, blank=True)
    has_publication = models.BooleanField(default=False)
    publication_title = models.TextField(default="")
    publication_doi = models.CharField(max_length=64, blank=True)
    publication_authors = ArrayField(models.TextField(), default=[])
    pubmed_id = models.CharField(max_length=32, blank=True)
    source_first_published = models.DateTimeField(null=True)
    source_last_modified = models.DateTimeField(null=True)

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
        return super(Experiment, self).save(*args, **kwargs)

    def to_metadata_dict(self):
        """ Render this Experiment as a dict """

        metadata = {}
        metadata['title'] = self.title
        metadata['accession_code'] = self.accession_code
        metadata['organisms'] = [organism.name for organism in self.organisms.all()]
        metadata['description'] = self.description
        metadata['protocol_description'] = self.protocol_description
        metadata['technology'] = self.technology
        metadata['submitter_institution'] = self.submitter_institution
        metadata['has_publication'] = self.has_publication
        metadata['publication_title'] = self.publication_title
        metadata['publication_doi'] = self.publication_doi
        metadata['pubmed_id'] = self.pubmed_id
        if self.source_first_published:
            metadata['source_first_published'] = self.source_first_published.strftime(
                '%Y-%m-%dT%H:%M:%S')
        else:
            metadata['source_first_published'] = ''
        if self.source_last_modified:
            metadata['source_last_modified'] = self.source_last_modified.strftime(
                '%Y-%m-%dT%H:%M:%S')
        else:
            metadata['source_last_modified'] = ''

        return metadata

    @property
    def platforms(self):
        """ Returns a list of related pipelines """
        return [p for p in self.samples.values_list('platform_name', flat=True).distinct()]

    @property
    def pretty_platforms(self):
        """ Returns a prettified list of related pipelines """
        return list(set([p.pretty_platform for p in self.samples.exclude(platform_name__exact='')]))

    def get_processed_samples(self):
        return self.samples.filter(is_processed=True)

class ExperimentAnnotation(models.Model):
    """ Semi-standard information associated with an Experiment """

    class Meta:
        db_table = "experiment_annotations"
        base_manager_name = 'public_objects'

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=False)

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
        return super(ExperimentAnnotation, self).save(*args, **kwargs)

class Pipeline(models.Model):
    """Pipeline that is associated with a series of ComputationalResult records."""

    name = models.CharField(max_length=255)
    steps = ArrayField(models.IntegerField(), default=[])

    class Meta:
        db_table = "pipelines"


class Processor(models.Model):
    """Processor associated with a certain ComputationalResult."""

    name = models.CharField(max_length=255)
    version = models.CharField(max_length=64)
    docker_image = models.CharField(max_length=255)
    environment = JSONField(default={})

    class Meta:
        db_table = "processors"
        unique_together = ('name', 'version')

    def __str__(self):
        return "Processor: %s (version: %s, docker_image: %s)" % (self.name, self.version, self.docker_image)


class ComputationalResult(models.Model):
    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"
        base_manager_name = 'public_objects'

    def __str__(self):
        return "ComputationalResult " + str(self.pk) + ": " + str(self.pipeline)

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    commands = ArrayField(models.TextField(), default=[])
    processor = models.ForeignKey(Processor, blank=True, null=True, on_delete=models.CASCADE)
    is_ccdl = models.BooleanField(default=True)
    # TODO: "pipeline" field is now redundant due to "processor". Should be removed later.
    # Human-readable nickname for this computation
    pipeline = models.CharField(max_length=255)

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


class ComputationalResultAnnotation(models.Model):
    """ Non-standard information associated with an ComputationalResult """

    class Meta:
        db_table = "computational_result_annotations"
        base_manager_name = 'public_objects'

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Properties
    data = HStoreField(default={})
    is_ccdl = models.BooleanField(default=True)

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
        base_manager_name = 'public_objects'

    def __str__(self):
        return "OrganismIndex " + str(self.pk) + ": " + self.organism.name + ' [' + self.index_type + '] - ' + str(self.salmon_version)

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # ex., "TRANSCRIPTOME_LONG", "TRANSCRIPTOME_SHORT"
    index_type = models.CharField(max_length=255)

    # This corresponds to Ensembl's release number:
    # http://ensemblgenomes.org/info/about/release_cycle
    # Determined by hitting:
    # http://rest.ensembl.org/info/software?content-type=application/json
    source_version = models.CharField(max_length=255, default="93")

    # The name of the genome assembly used which corresponds to 'GRCh38' in:
    # ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    assembly_name = models.CharField(max_length=255, default="UNKNOWN")

    # This matters, for instance salmon 0.9.0 indexes don't work with 0.10.0
    salmon_version = models.CharField(max_length=255, default="0.9.1")

    # S3 Information
    s3_url = models.CharField(max_length=255, default="")

    # Common Properties
    is_public = models.BooleanField(default=True)
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def upload_to_s3(self, absolute_file_path, bucket_name, logger):
        if bucket_name is not None:
            s3_key = self.organism.name + '_' + self.index_type + '.tar.gz'
            S3.upload_file(absolute_file_path, bucket_name, s3_key,
                           ExtraArgs={'ACL': 'public-read'})
            self.s3_url = get_s3_url(bucket_name, s3_key)
            logger.info("Upload complete")
        else:
            logger.info("No S3 bucket in environment, not uploading")

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

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # File Properties
    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField(blank=True, null=True)
    sha1 = models.CharField(max_length=64)

    # AWS
    s3_bucket = models.CharField(max_length=255, blank=True, null=True)
    s3_key = models.CharField(max_length=255, blank=True, null=True)

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

    def sync_to_s3(self, s3_bucket=None, s3_key=None) -> bool:
        """ Syncs this OriginalFile to AWS S3.
        """
        if not settings.RUNNING_IN_CLOUD:
            return True

        self.s3_bucket = s3_bucket
        self.s3_key = s3_key

        try:
            S3.upload_file(
                        self.absolute_file_path,
                        s3_bucket,
                        s3_key,
                        ExtraArgs={
                            'ACL': 'public-read',
                            'StorageClass': 'STANDARD_IA'
                        }
                    )
            self.save()
        except Exception as e:
            logger.exception(e, original_file_id=self.pk)
            return False

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

    def get_display_name(self):
        """ For dev convenience """
        if not self.filename:
            return self.source_filename
        else:
            return self.filename

    def sync_from_s3(self):
        """ Downloads a file from S3 to the local file system.
        Returns the absolute file path.
        """
        if not settings.RUNNING_IN_CLOUD:
            return self.absolute_file_path

        try:
            S3.download_file(
                        self.s3_bucket,
                        self.s3_key,
                        self.absolute_file_path
                    )
            return self.absolute_file_path
        except Exception as e:
            logger.exception(e, original_file_id=self.pk)
            return None

    def get_synced_file_path(self):
        """ Fetches the absolute file path to this ComputedFile, fetching from S3 if it
        isn't already available locally. """
        if os.path.exists(self.absolute_file_path):
            return self.absolute_file_path
        else:
            return self.sync_from_s3()

    def delete_local_file(self):
        """ Deletes this file from the local file system."""
        if not settings.RUNNING_IN_CLOUD:
            return

        try:
            os.remove(self.absolute_file_path)
        except OSError:
            pass
        self.is_downloaded = False
        self.save()


class ComputedFile(models.Model):
    """ A representation of a file created by a data-refinery process """

    class Meta:
        db_table = "computed_files"

    def __str__(self):
        return "ComputedFile: " + str(self.filename)

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # File related
    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    # Relations
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Scientific
    is_smashable = models.BooleanField(default=False)
    is_qc = models.BooleanField(default=False)

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

    def sync_to_s3(self, s3_bucket=None, s3_key=None) -> bool:
        """ Syncs a file to AWS S3.
        """
        if not settings.RUNNING_IN_CLOUD:
            return True

        self.s3_bucket = s3_bucket
        self.s3_key = s3_key

        try:
            S3.upload_file(
                        self.absolute_file_path,
                        s3_bucket,
                        s3_key,
                        ExtraArgs={
                            'ACL': 'public-read',
                            'StorageClass': 'STANDARD_IA'
                        }
                    )
            self.save()
        except Exception as e:
            logger.exception(e, computed_file_id=self.pk)
            return False

        return True

    def sync_from_s3(self, force=False):
        """ Downloads a file from S3 to the local file system.
        Returns the absolute file path.
        """
        if not settings.RUNNING_IN_CLOUD and not force:
            return self.absolute_file_path

        try:
            S3.download_file(
                        self.s3_bucket,
                        self.s3_key,
                        self.absolute_file_path
                    )
            return self.absolute_file_path
        except Exception as e:
            logger.exception(e, computed_file_id=self.pk)
            return None

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

    def delete_local_file(self):
        """ Deletes a file from the path and actually removes it from the file system,
        resetting the is_downloaded flag to false. Can be refetched from source if needed. """
        if not settings.RUNNING_IN_CLOUD:
            return

        try:
            os.remove(self.absolute_file_path)
        except OSError:
            pass

    def get_synced_file_path(self, force=False):
        """ Fetches the absolute file path to this ComputedFile, fetching from S3 if it
        isn't already available locally. """
        if os.path.exists(self.absolute_file_path):
            return self.absolute_file_path
        else:
            return self.sync_from_s3(force)

class Dataset(models.Model):
    """ A Dataset is a desired set of experiments/samples to smash and download """

    AGGREGATE_CHOICES = (
        ('ALL', 'All'),
        ('EXPERIMENT', 'Experiment'),
        ('SPECIES', 'Species')
    )

    SCALE_CHOICES = (
        ('NONE', 'None'),
        ('MINMAX', 'Minmax'),
        ('STANDARD', 'Standard'),
        ('ROBUST', 'Robust'),
    )

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Experiments and samples live here: {'E-ABC-1': ['SAMP1', 'SAMP2']}
    # This isn't going to be queryable, so we can use JSON-in-text, just make
    # sure we validate properly in and out!
    data = JSONField(default={})

    # Processing properties
    aggregate_by = models.CharField(max_length=255, choices=AGGREGATE_CHOICES, default="EXPERIMENT")
    scale_by = models.CharField(max_length=255, choices=SCALE_CHOICES, default="MINMAX")

    # State properties
    is_processing = models.BooleanField(default=False)  # Data is still editable when False
    is_processed = models.BooleanField(default=False)  # Result has been made
    is_available = models.BooleanField(default=False)  # Result is ready for delivery

    # Fail handling
    success = models.NullBooleanField(null=True)
    failure_reason = models.TextField()

    # Delivery properties
    email_address = models.CharField(max_length=255, blank=True, null=True)
    email_sent = models.BooleanField(default=False)  # Result has been made
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

    def get_samples(self):
        """ Retuns all of the Sample objects in this Dataset """

        all_samples = []
        for sample_list in self.data.values():
            all_samples = all_samples + sample_list
        all_samples = list(set(all_samples))

        return Sample.objects.filter(accession_code__in=all_samples)

    def get_experiments(self):
        """ Retuns all of the Experiments objects in this Dataset """

        all_experiments = []
        for experiment in self.data.keys():
            all_experiments.append(experiment)
        all_experiments = list(set(all_experiments))

        return Experiment.objects.filter(accession_code__in=all_experiments)

    def get_samples_by_experiment(self):
        """ Returns a dict of sample QuerySets, for samples grouped by experiment. """
        all_samples = {}

        for experiment, samples in self.data.items():
            all_samples[experiment] = Sample.objects.filter(accession_code__in=samples)

        return all_samples

    def get_samples_by_species(self):
        """ Returns a dict of sample QuerySets, for samples grouped by species. """

        by_species = {}
        all_samples = self.get_samples()
        for sample in all_samples:
            if not by_species.get(sample.organism.name, None):
                by_species[sample.organism.name] = [sample]
            else:
                by_species[sample.organism.name].append(sample)

        return by_species

    def get_aggregated_samples(self):
        """ Uses aggregate_by to return smasher-ready a sample dict. """

        if self.aggregate_by == "ALL":
            return {'ALL': self.get_samples()}
        elif self.aggregate_by == "EXPERIMENT":
            return self.get_samples_by_experiment()
        else:
            return self.get_samples_by_species()

    def is_cross_technology(self):
        """ Determine if this involves both Microarray + RNASeq"""

        if len(self.get_samples().values('technology').distinct()) > 1:
            return True
        else:
            return False

    def s3_url(self):
        """ Render the resulting S3 URL """
        if (self.s3_key) and (self.s3_bucket):
            return "https://s3.amazonaws.com/" + self.s3_bucket + "/" + self.s3_key
        else:
            return None

class APIToken(models.Model):
    """ Required for starting a smash job """

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Activation
    is_activated = models.BooleanField(default=False)

    # Common Properties
    created_at = models.DateTimeField(editable=False, default=timezone.now)
    last_modified = models.DateTimeField(default=timezone.now)

    def save(self, *args, **kwargs):
        """ On save, update timestamps """
        current_time = timezone.now()
        if not self.id:
            self.created_at = current_time
        self.last_modified = current_time
        return super(APIToken, self).save(*args, **kwargs)

    @property
    def terms_and_conditions(self):
        """ """
        return settings.TERMS_AND_CONDITIONS

"""
# Associations

These represent the relationships between items in the other tables.
"""


class ExperimentSampleAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_sample_associations"
        unique_together = ('experiment', 'sample')


class ExperimentOrganismAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_organism_associations"
        unique_together = ('experiment', 'organism')


class DownloaderJobOriginalFileAssociation(models.Model):

    downloader_job = models.ForeignKey(
        "data_refinery_common.DownloaderJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "downloaderjob_originalfile_associations"
        unique_together = ('downloader_job', 'original_file')


class ProcessorJobOriginalFileAssociation(models.Model):

    processor_job = models.ForeignKey(
        "data_refinery_common.ProcessorJob", blank=False, null=False, on_delete=models.CASCADE)
    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "processorjob_originalfile_associations"
        unique_together = ('processor_job', 'original_file')


class ProcessorJobDatasetAssociation(models.Model):

    processor_job = models.ForeignKey(
        "data_refinery_common.ProcessorJob", blank=False, null=False, on_delete=models.CASCADE)
    dataset = models.ForeignKey(Dataset, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "processorjob_dataset_associations"


class OriginalFileSampleAssociation(models.Model):

    original_file = models.ForeignKey(
        OriginalFile, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "original_file_sample_associations"
        unique_together = ('original_file', 'sample')


class SampleResultAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "sample_result_associations"
        unique_together = ('result', 'sample')


class SampleComputedFileAssociation(models.Model):

    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)
    computed_file = models.ForeignKey(
        ComputedFile, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "sample_computed_file_associations"
        unique_together = ('sample', 'computed_file')

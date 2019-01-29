import hashlib
import io
import os
import shutil
import pytz
import uuid
import boto3

from botocore.client import Config

from datetime import datetime
from functools import partial
from typing import Dict, Set

from django.conf import settings
from django.contrib.postgres.fields import ArrayField, JSONField
from django.db import transaction
from django.db import models
from django.utils import timezone

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.organism import Organism
from data_refinery_common.utils import get_env_variable, get_s3_url, calculate_file_size, calculate_sha1

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
        base_manager_name = "public_objects"
        get_latest_by = "created_at"

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
    experiments = models.ManyToManyField('Experiment', through='ExperimentSampleAssociation')

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
    protocol_info = JSONField(default=dict)

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
        """Render this Sample as a dict."""
        metadata = {}
        metadata['refinebio_title'] = self.title
        metadata['refinebio_accession_code'] = self.accession_code
        metadata['refinebio_organism'] = self.organism.name if self.organism else None
        metadata['refinebio_source_database'] = self.source_database
        metadata['refinebio_source_archive_url'] = self.source_archive_url
        metadata['refinebio_sex'] = self.sex
        metadata['refinebio_age'] = self.age or ''
        metadata['refinebio_specimen_part'] = self.specimen_part
        metadata['refinebio_genetic_information'] = self.genotype
        metadata['refinebio_disease'] = self.disease
        metadata['refinebio_disease_stage'] = self.disease_stage
        metadata['refinebio_cell_line'] = self.cell_line
        metadata['refinebio_treatment'] = self.treatment
        metadata['refinebio_race'] = self.race
        metadata['refinebio_subject'] = self.subject
        metadata['refinebio_compound'] = self.compound
        metadata['refinebio_time'] = self.time
        metadata['refinebio_platform'] = self.pretty_platform
        metadata['refinebio_annotations'] = [
            data for data in self.sampleannotation_set.all().values_list('data', flat=True)
        ]

        return metadata

    # Returns a set of ProcessorJob objects but we cannot specify
    # that in type hints because it hasn't been declared yet.
    def get_processor_jobs(self) -> Set:
        processor_jobs = set()
        for original_file in self.original_files.all():
            for processor_job in original_file.processor_jobs.all():
                processor_jobs.add(processor_job)

        return processor_jobs

    def get_result_files(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return self.computed_files.all()

    def get_most_recent_smashable_result_file(self):
        """ Get the most recent of the ComputedFile objects associated with this Sample """
        try:
            return self.computed_files.filter(
                            is_smashable=True
                        ).latest()
        except Exception as e:
            # This sample has no smashable files yet.
            return None

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
    data = JSONField(default=dict)
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
    alternate_accession_code = models.CharField(max_length=64, unique=True, null=True)

    # Historical Properties
    source_database = models.CharField(max_length=32)  # "ArrayExpress, "SRA", "GEO"
    source_url = models.TextField()

    # Properties
    # I was always under the impression that TextFields were slower
    # than CharFields, however the Postgres documentation disagrees:
    # https://www.postgresql.org/docs/9.0/static/datatype-character.html
    title = models.TextField()
    description = models.TextField()
    protocol_description = JSONField(default=dict)
    technology = models.CharField(max_length=256, blank=True)
    submitter_institution = models.CharField(max_length=256, blank=True)
    has_publication = models.BooleanField(default=False)
    publication_title = models.TextField(default="")
    publication_doi = models.CharField(max_length=64, blank=True)
    publication_authors = ArrayField(models.TextField(), default=list)
    pubmed_id = models.CharField(max_length=32, blank=True)
    source_first_published = models.DateTimeField(null=True)
    source_last_modified = models.DateTimeField(null=True)

    # Cached Computed Properties
    num_total_samples = models.IntegerField(default=0)
    num_processed_samples = models.IntegerField(default=0)
    sample_metadata_fields = ArrayField(models.TextField(), default=list)
    organism_names = ArrayField(models.TextField(), default=list)
    platform_names = ArrayField(models.TextField(), default=list)

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

        if self.accession_code and not self.alternate_accession_code:
            if self.accession_code.startswith('GSE'):
                self.alternate_accession_code = 'E-GEOD-' + self.accession_code[3:]
            elif self.accession_code.startswith('E-GEOD-'):
                self.alternate_accession_code = 'GSE' + self.accession_code[7:]

        return super(Experiment, self).save(*args, **kwargs)

    def update_num_samples(self):
        """ Update our cache values """
        self.num_total_samples = self.samples.count()
        self.num_processed_samples = self.samples.filter(is_processed=True).count()
        self.save()

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

    def get_sample_metadata_fields(self):
        """ Get all metadata fields that are non-empty for at least one sample in the experiment.
        See https://github.com/AlexsLemonade/refinebio-frontend/issues/211 for why this is needed.
        """
        fields = []

        possible_fields = ['sex', 'age', 'specimen_part', 'genotype', 'disease', 'disease_stage',
                           'cell_line', 'treatment', 'race', 'subject', 'compound', 'time']
        samples = self.samples.all()
        for field in possible_fields:
            for sample in samples:
                if getattr(sample, field) != None and getattr(sample, field) != '':
                    fields.append(field)
                    break

        return fields

    def update_sample_metadata_fields(self):
        self.sample_metadata_fields = self.get_sample_metadata_fields()

    def update_organism_names(self):
        self.organism_names = self.get_organism_names()

    def update_platform_names(self):
        self.platform_names = self.get_platform_names()

    def get_sample_technologies(self):
        """ Get a list of unique technologies for all of the associated samples
        """
        return list(set([sample.technology for sample in self.samples.all()]))

    def get_organism_names(self):
        """ Get a list of unique technologies for all of the associated samples
        """
        return list(set([organism.name for organism in self.organisms.all()]))

    def get_platform_names(self):
        """ Get a list of unique platforms for all of the associated samples
        """
        return list(set([sample.platform_name for sample in self.samples.all()]))

    @property
    def platforms(self):
        """ Returns a list of related pipelines """
        return list(set([sample.platform_name for sample in self.samples.all()]))

    @property
    def pretty_platforms(self):
        """ Returns a prettified list of related pipelines """
        return list(set([sample.pretty_platform for sample in self.samples.all()]))

    @property
    def processed_samples(self):
        return list([sample.accession_code for sample in self.samples.all() if sample.is_processed == True])

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
    data = JSONField(default=dict)
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
    steps = ArrayField(models.IntegerField(), default=list)

    class Meta:
        db_table = "pipelines"


class Processor(models.Model):
    """Processor associated with a certain ComputationalResult."""

    name = models.CharField(max_length=255)
    version = models.CharField(max_length=64)
    docker_image = models.CharField(max_length=255)
    environment = JSONField(default=dict)

    class Meta:
        db_table = "processors"
        unique_together = ('name', 'version', 'docker_image', 'environment')

    def __str__(self):
        return "Processor: %s (version: %s, docker_image: %s)" % (self.name, self.version, self.docker_image)


class ComputationalResult(models.Model):
    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"
        base_manager_name = 'public_objects'

    def __str__(self):
        processor_name_str = ""
        if self.processor:
            processor_name_str = ": " + str(self.processor.name)

        return "ComputationalResult " + str(self.pk) + processor_name_str

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    commands = ArrayField(models.TextField(), default=list)
    processor = models.ForeignKey(Processor, blank=True, null=True, on_delete=models.CASCADE)

    samples = models.ManyToManyField('Sample', through='SampleResultAssociation')

    # The Organism Index used to process the sample.
    organism_index = models.ForeignKey('OrganismIndex', blank=True, null=True, on_delete=models.SET_NULL)

    is_ccdl = models.BooleanField(default=True)

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
    data = JSONField(default=dict)
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

    # We keep the director unextracted on the shared filesystem so all
    # Salmon jobs can access it.
    absolute_directory_path = models.CharField(max_length=255, blank=True, null=True, default="")

    # S3 Information
    s3_url = models.CharField(max_length=255, default="")

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
    processor_jobs = models.ManyToManyField('data_refinery_common.ProcessorJob', through='ProcessorJobOriginalFileAssociation')
    downloader_jobs = models.ManyToManyField('data_refinery_common.DownloaderJob', through='DownloaderJobOriginalFileAssociation')

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

    def get_display_name(self):
        """ For dev convenience """
        if not self.filename:
            return self.source_filename
        else:
            return self.filename

    def delete_local_file(self):
        """ Deletes this file from the local file system."""
        try:
            os.remove(self.absolute_file_path)
        except OSError:
            pass
        except TypeError:
            pass
        except Exception as e:
            logger.exception("Unexpected delete file exception.",
                absolute_file_path=self.absolute_file_path
            )
        self.is_downloaded = False
        self.save()


class ComputedFile(models.Model):
    """ A representation of a file created by a data-refinery process """

    class Meta:
        db_table = "computed_files"
        get_latest_by = "created_at"

    def __str__(self):
        return "ComputedFile: " + str(self.filename)

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Object relations
    samples = models.ManyToManyField('Sample', through='SampleComputedFileAssociation')

    # File related
    filename = models.CharField(max_length=255)
    absolute_file_path = models.CharField(max_length=255, blank=True, null=True)
    # TODO: make this work w/ migrations:
    # absolute_file_path = models.CharField(max_length=255)
    size_in_bytes = models.BigIntegerField()
    sha1 = models.CharField(max_length=64)

    # Relations
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    # Scientific
    is_smashable = models.BooleanField(default=False)
    is_qc = models.BooleanField(default=False)
    is_qn_target = models.BooleanField(default=False)

    # Compendia details
    is_compendia = models.BooleanField(default=False)
    compendia_organism = models.ForeignKey(Organism,
                                        blank=True,
                                        null=True, 
                                        on_delete=models.CASCADE
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

    def sync_to_s3(self, s3_bucket=None, s3_key=None) -> bool:
        """ Syncs a file to AWS S3.
        """
        if not settings.RUNNING_IN_CLOUD:
            return False

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

        target_directory = os.path.dirname(self.absolute_file_path)
        os.makedirs(target_directory, exist_ok=True)

        try:
            S3.download_file(
                        self.s3_bucket,
                        self.s3_key,
                        path
                    )

            # Veryify sync integrity
            synced_sha1 = calculate_sha1(path)

            if self.sha1 != synced_sha1:
                raise AssertionError("SHA1 of downloaded ComputedFile doesn't match database SHA1!")

            return path
        except Exception as e:
            logger.exception(e, computed_file_id=self.pk)
            return None

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
        except Exception as e:
            logger.exception("Unexpected delete file exception.",
                absolute_file_path=self.absolute_file_path
            )

    def delete_s3_file(self, force=False):
        # If we're not running in the cloud then we shouldn't try to
        # delete something from S3 unless force is set.
        if not settings.RUNNING_IN_CLOUD and not force:
            return False

        try:
            S3.delete_object(Bucket=self.s3_bucket, Key=self.s3_key)
            return True
        except:
            logger.exception("Failed to delete S3 object for Computed File.",
                             computed_file=self.id,
                             s3_object=self.s3_key)
            return False

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
        """ Render the resulting S3 URL """
        return self.get_s3_url()

    def get_s3_url(self):
        """ Render the resulting S3 URL """
        if (self.s3_key) and (self.s3_bucket):
            return "https://s3.amazonaws.com/" + self.s3_bucket + "/" + self.s3_key
        else:
            return None

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
    data = JSONField(default=dict)

    # Processing properties
    aggregate_by = models.CharField(max_length=255, choices=AGGREGATE_CHOICES, default="EXPERIMENT")
    scale_by = models.CharField(max_length=255, choices=SCALE_CHOICES, default="NONE")
    quantile_normalize = models.BooleanField(default=True)

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

    size_in_bytes = models.BigIntegerField(blank=True, null=True, default=0)
    sha1 = models.CharField(max_length=64, null=True, default='')

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
        all_experiments = self.data.keys()
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
        """ Uses aggregate_by to return a smasher-ready sample dict. """

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

    @property
    def has_email(self):
        """ Returns if the email is set or not """
        return bool(self.email_address)

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


class ExperimentResultAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_result_associations"
        unique_together = ('result', 'experiment')

import hashlib
import io
import os
import pytz
import uuid

from datetime import datetime
from functools import partial

from django.contrib.postgres.fields import ArrayField, HStoreField, JSONField
from django.db import transaction
from django.db import models
from django.utils import timezone

from data_refinery_common.models.organism import Organism

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
        return metadata

    def get_result_files(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return ComputedFile.objects.filter(result__in=self.results.all())

    def get_most_recent_smashable_result_file(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return ComputedFile.objects.filter(
                        result__in=self.results.all(),
                        is_smashable=True,
                    ).first()

    @property
    def pipelines(self):
        """ Returns a list of related pipelines """
        return [p for p in self.results.values_list('pipeline', flat=True).distinct()]

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


class ComputationalResult(models.Model):
    """ Meta-information about the output of a computer process. (Ex Salmon) """

    class Meta:
        db_table = "computational_results"
        base_manager_name = 'public_objects'

    def __str__(self):
        return "ComputationalResult: " + str(self.pk)

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    commands = ArrayField(models.TextField())
    is_ccdl = models.BooleanField(default=True)

    # TODO: This field should be changed into "processor" (foreign key to "Processor").
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


class Pipeline(models.Model):
    "Pipeline that is associated with a series of ComputationalResult records."""

    name = models.CharField(max_length=255)
    steps = ArrayField(models.IntegerField(), default=[])

    class Meta:
        db_table = "pipelines"


class Processor(models.Model):
    """Processor associated with a certain ComputationalResult."""

    name = models.CharField(max_length=255)
    docker_image = models.CharField(max_length=255)
    environment = JSONField(default={})

    class Meta:
        db_table = "processors"
        unique_together = ('name', 'docker_image')

    def __str__(self):
        return "Processor: %s (docker_image: %s" % (name, docker_image)

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

    # Managers
    objects = models.Manager()
    public_objects = PublicObjectsManager()

    # Relations
    organism = models.ForeignKey(Organism, blank=False, null=False, on_delete=models.CASCADE)

    # ex., "TRANSCRIPTOME_LONG", "TRANSCRIPTOME_SHORT"
    index_type = models.CharField(max_length=255)
    # This corresponds to Ensembl's release number:
    # http://ensemblgenomes.org/info/about/release_cycle
    # Determined by hitting:
    # http://rest.ensembl.org/info/software?content-type=application/json
    source_version = models.CharField(max_length=255)
    result = models.ForeignKey(
        ComputationalResult, blank=False, null=False, on_delete=models.CASCADE)

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
    s3_bucket = models.CharField(max_length=255)
    s3_key = models.CharField(max_length=255)

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

from typing import Set

from django.db import models
from django.utils import timezone

from data_refinery_common.models.computed_file import ComputedFile
from data_refinery_common.models.managers import ProcessedObjectsManager, PublicObjectsManager


class Sample(models.Model):
    """
    An individual sample.
    """

    class Meta:
        db_table = "samples"
        base_manager_name = "public_objects"
        get_latest_by = "created_at"

        indexes = [
            models.Index(fields=["accession_code"]),
        ]

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
    organism = models.ForeignKey("Organism", blank=True, null=True, on_delete=models.SET_NULL)
    results = models.ManyToManyField("ComputationalResult", through="SampleResultAssociation")
    original_files = models.ManyToManyField("OriginalFile", through="OriginalFileSampleAssociation")
    computed_files = models.ManyToManyField("ComputedFile", through="SampleComputedFileAssociation")
    experiments = models.ManyToManyField("Experiment", through="ExperimentSampleAssociation")

    # Historical Properties
    source_database = models.CharField(max_length=255, blank=False)
    source_archive_url = models.CharField(max_length=255)
    source_filename = models.CharField(max_length=255, blank=False)
    source_absolute_file_path = models.CharField(max_length=255)
    has_raw = models.BooleanField(default=True)  # Did this sample have a raw data source?

    # Technological Properties
    platform_accession_code = models.CharField(max_length=256, blank=True)
    platform_name = models.CharField(max_length=256, blank=True)
    technology = models.CharField(max_length=256, blank=True)  # MICROARRAY, RNA-SEQ
    manufacturer = models.CharField(max_length=256, blank=True)
    protocol_info = models.JSONField(default=dict)

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

    # Blacklisting
    is_blacklisted = models.BooleanField(default=False)

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

    def to_metadata_dict(self, computed_file=None):
        """Render this Sample as a dict."""
        metadata = {}
        metadata["refinebio_title"] = self.title
        metadata["refinebio_accession_code"] = self.accession_code
        metadata["refinebio_organism"] = self.organism.name if self.organism else None
        metadata["refinebio_source_database"] = self.source_database
        metadata["refinebio_source_archive_url"] = self.source_archive_url
        metadata["refinebio_sex"] = self.sex
        metadata["refinebio_age"] = self.age or ""
        metadata["refinebio_specimen_part"] = self.specimen_part
        metadata["refinebio_genetic_information"] = self.genotype
        metadata["refinebio_disease"] = self.disease
        metadata["refinebio_disease_stage"] = self.disease_stage
        metadata["refinebio_cell_line"] = self.cell_line
        metadata["refinebio_treatment"] = self.treatment
        metadata["refinebio_race"] = self.race
        metadata["refinebio_subject"] = self.subject
        metadata["refinebio_compound"] = self.compound
        metadata["refinebio_time"] = self.time
        metadata["refinebio_platform"] = self.pretty_platform
        metadata["refinebio_processed"] = self.has_raw
        metadata["refinebio_annotations"] = [
            data for data in self.sampleannotation_set.all().values_list("data", flat=True)
        ]

        if computed_file and computed_file.result and computed_file.result.processor:
            metadata["refinebio_processor_id"] = computed_file.result.processor.id
            metadata["refinebio_processor_name"] = computed_file.result.processor.name
            metadata["refinebio_processor_version"] = computed_file.result.processor.version

        if self.attributes.count() > 0:
            metadata["other_metadata"] = [
                attribute.to_dict() for attribute in self.attributes.all()
            ]

        return metadata

    # Returns a set of ProcessorJob objects but we cannot specify
    # that in type hints because it hasn't been declared yet.
    def get_processor_jobs(self) -> Set:
        processor_jobs = set()
        for original_file in self.original_files.prefetch_related("processor_jobs").all():
            for processor_job in original_file.processor_jobs.all():
                processor_jobs.add(processor_job)

        return processor_jobs

    # Returns a set of DownloaderJob objects but we cannot specify
    # that in type hints because it hasn't been declared yet.
    def get_downloader_jobs(self) -> Set:
        downloader_jobs = set()
        for original_file in self.original_files.prefetch_related("downloader_jobs").all():
            for downloader_job in original_file.downloader_jobs.all():
                downloader_jobs.add(downloader_job)

        return downloader_jobs

    def get_result_files(self):
        """ Get all of the ComputedFile objects associated with this Sample """
        return self.computed_files.all()

    def get_most_recent_smashable_result_file(self):
        """ Get the most recent of the ComputedFile objects associated with this Sample """
        try:
            latest_computed_file = self.computed_files.filter(
                is_public=True, is_smashable=True, s3_bucket__isnull=False, s3_key__isnull=False
            ).latest()
            return latest_computed_file
        except ComputedFile.DoesNotExist as e:
            # This sample has no smashable files yet.
            return None

    def get_most_recent_quant_sf_file(self):
        """ Returns the latest quant.sf file that was generated for this sample.
        Note: We don't associate that file to the computed_files of this sample, that's
        why we have to go through the computational results. """
        return (
            ComputedFile.objects.filter(
                result__in=self.results.all(),
                filename="quant.sf",
                s3_key__isnull=False,
                s3_bucket__isnull=False,
            )
            .order_by("-created_at")
            .first()
        )

    @property
    def pretty_platform(self):
        """ Turns

        [HT_HG-U133_Plus_PM] Affymetrix HT HG-U133+ PM Array Plate

        into

        Affymetrix HT HG-U133+ PM Array Plate (hthgu133pluspm)

        """
        if "]" in self.platform_name:
            platform_base = self.platform_name.split("]")[1].strip()
        else:
            platform_base = self.platform_name
        return platform_base + " (" + self.platform_accession_code + ")"

    @property
    def experiment_accession_codes(self):
        return [e.accession_code for e in self.experiments.all()]

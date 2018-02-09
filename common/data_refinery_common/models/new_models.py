from django.db import transaction
from django.db import models
from django.utils import timezone
from datetime import datetime

import pytz

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
    source_filename = models.CharField(max_length=255)

    # Scientific Properties
    has_raw = models.BooleanField(default=True) # Did this sample have a raw data source?
    has_derived = models.BooleanField(default=False) # Did this sample have a pre-derived data source?

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

class ComputationalResult(models.Model):

    class Meta:
        db_table = "computational_results"

class CompultationalResultAnnotation(models.Model):

    class Meta:
        db_table = "computational_result_annotations"

class Gene(models.Model):

    class Meta:
        db_table = "genes"

"""
# Associations

These represent the relationships between items in the other tables. 
"""

class ExperimentSampleAssociation(models.Model):

    experiment = models.ForeignKey(Experiment, blank=False, null=False, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = "experiment_sample_associations"

class SampleResultAssociation(models.Model):

    class Meta:
        db_table = "sample_result_associations"

class ExperimentResultAssociation(models.Model):

    class Meta:
        db_table = "experiment_result_associations"

"""
# Files

These are the database representations of files
which live on local disk, on ephemeral storage,
or on AWS cloud services.
"""

"""
# Jobs

These models correspond to the status of various task workers
responsible for the downloading and transformation of our
input data.

"""


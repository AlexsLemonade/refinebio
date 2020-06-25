import uuid

from django.conf import settings
from django.contrib.postgres.fields import JSONField
from django.db import models
from django.utils import timezone

import boto3
from botocore.client import Config

from data_refinery_common.models.experiment import Experiment
from data_refinery_common.models.sample import Sample

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))


class Dataset(models.Model):
    """ A Dataset is a desired set of experiments/samples to smash and download """

    AGGREGATE_CHOICES = (("ALL", "All"), ("EXPERIMENT", "Experiment"), ("SPECIES", "Species"))

    SCALE_CHOICES = (
        ("NONE", "None"),
        ("MINMAX", "Minmax"),
        ("STANDARD", "Standard"),
        ("ROBUST", "Robust"),
    )

    SVD_ALGORITHM_CHOICES = (
        ("NONE", "None"),
        ("RANDOMIZED", "randomized"),
        ("ARPACK", "arpack"),
    )

    # ID
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)

    # Experiments and samples live here: {'E-ABC-1': ['SAMP1', 'SAMP2']}
    # This isn't going to be queryable, so we can use JSON-in-text, just make
    # sure we validate properly in and out!
    data = JSONField(
        default=dict,
        help_text=(
            "This is a dictionary where the keys are experiment accession codes and the values are"
            " lists with sample accession codes. Eg: `{'E-ABC-1': ['SAMP1', 'SAMP2']}`"
        ),
    )

    # Processing properties
    aggregate_by = models.CharField(
        max_length=255,
        choices=AGGREGATE_CHOICES,
        default="EXPERIMENT",
        help_text=(
            "Specifies how samples are"
            " [aggregated](http://docs.refine.bio/en/latest/main_text.html#aggregations)."
        ),
    )
    scale_by = models.CharField(
        max_length=255,
        choices=SCALE_CHOICES,
        default="NONE",
        help_text=(
            "Specifies options for"
            " [transformations](http://docs.refine.bio/en/latest/main_text.html#transformations)."
        ),
    )
    quantile_normalize = models.BooleanField(
        default=True,
        help_text=(
            "Part of the advanced options. Allows "
            "[skipping quantile normalization](http://docs.refine.bio/en/latest/faq.html#what-does-it-mean-to-skip-quantile-normalization-for-rna-seq-samples)"  # noqa
            " for RNA-Seq samples."
        ),
    )
    quant_sf_only = models.BooleanField(
        default=False, help_text="Include only quant.sf files in the generated dataset."
    )
    svd_algorithm = models.CharField(
        max_length=255,
        choices=SVD_ALGORITHM_CHOICES,
        default="NONE",
        help_text="Specifies choice of SVD algorithm",
    )

    # State properties
    is_processing = models.BooleanField(default=False)  # Data is still editable when False
    is_processed = models.BooleanField(default=False)  # Result has been made
    is_available = models.BooleanField(default=False)  # Result is ready for delivery

    processor_jobs = models.ManyToManyField(
        "data_refinery_common.ProcessorJob", through="ProcessorJobDataSetAssociation"
    )

    # Fail handling
    success = models.NullBooleanField(null=True)
    failure_reason = models.TextField()

    # Delivery properties
    email_address = models.EmailField(max_length=255, blank=True, null=True)
    email_ccdl_ok = models.BooleanField(default=False)
    expires_on = models.DateTimeField(blank=True, null=True)

    # Deliverables
    s3_bucket = models.CharField(max_length=255)
    s3_key = models.CharField(max_length=255)

    size_in_bytes = models.BigIntegerField(
        blank=True,
        null=True,
        default=0,
        help_text="Contains the size in bytes of the processed dataset.",
    )
    sha1 = models.CharField(max_length=64, null=True, default="")

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

    def get_total_samples(self):
        """ Returns the total number of samples, this counts the number of unique
        accession codes in `data`. """
        return len(
            set(
                [
                    accession_code
                    for experiment in self.data.values()
                    for accession_code in experiment
                ]
            )
        )

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
            return {"ALL": self.get_samples()}
        elif self.aggregate_by == "EXPERIMENT":
            return self.get_samples_by_experiment()
        else:
            return self.get_samples_by_species()

    def is_cross_technology(self):
        """ Determine if this involves both Microarray + RNASeq"""

        if len(self.get_samples().values("technology").distinct()) > 1:
            return True
        else:
            return False

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

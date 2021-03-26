from typing import Dict

from django.test import TestCase
from django.utils import timezone

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentResultAssociation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    OriginalFile,
    OriginalFileSampleAssociation,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleResultAssociation,
)
from data_refinery_foreman.foreman.management.commands.retry_samples import retry_by_regex


def setup_experiment() -> Dict:
    """ Create an experiment with two samples where one of them has a processor job that failed
    because Batch restarted it."""

    # Create the experiment
    experiment_accession = "SRP095529"
    experiment = Experiment.objects.create(
        accession_code=experiment_accession, technology="RNA-SEQ"
    )

    zebrafish = Organism.get_object_for_name("DANIO_RERIO")

    accession_code = "S001"
    sample = Sample.objects.create(
        accession_code=accession_code,
        organism=zebrafish,
        source_database="SRA",
        technology="RNA-SEQ",
        platform_accession_code="IlluminaHiSeq1000",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    original_file = OriginalFile()
    original_file.filename = accession_code + ".SRA"
    original_file.source_filename = accession_code + ".SRA"
    original_file.save()

    OriginalFileSampleAssociation.objects.get_or_create(original_file=original_file, sample=sample)

    accession_code = "S002"
    sample = Sample.objects.create(
        accession_code=accession_code,
        organism=zebrafish,
        source_database="SRA",
        technology="RNA-SEQ",
        platform_accession_code="IlluminaHiSeq1000",
    )
    ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

    original_file = OriginalFile()
    original_file.filename = accession_code + ".SRA"
    original_file.source_filename = accession_code + ".SRA"
    original_file.save()

    OriginalFileSampleAssociation.objects.get_or_create(original_file=original_file, sample=sample)

    # add a failed processor job for the second sample
    processor_job = ProcessorJob()
    processor_job.start_time = timezone.now()
    processor_job.end_time = timezone.now()
    processor_job.no_retry = True
    processor_job.success = False
    processor_job.failure_reason = (
        "ProcessorJob has already completed with a fail - why are we here again?"
    )
    processor_job.save()

    processor_assoc = ProcessorJobOriginalFileAssociation()
    processor_assoc.original_file = original_file
    processor_assoc.processor_job = processor_job
    processor_assoc.save()

    return experiment


class RetrySamples(TestCase):
    """
    """

    def test(self):
        setup_experiment()
        retry_by_regex("ProcessorJob has already completed .*")

        dl_jobs = DownloaderJob.objects.all()
        self.assertEqual(dl_jobs.count(), 1)

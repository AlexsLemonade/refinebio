import datetime
import json
import os

from django.test import TestCase
from django.utils import timezone
from typing import Dict, List
from unittest.mock import MagicMock, Mock, patch, call

from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
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
    OriginalFile,
    OriginalFileSampleAssociation,
    Organism,
    OrganismIndex,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleResultAssociation,
)
from data_refinery_foreman.foreman.management.commands.rerun_salmon_old_samples import update_salmon_all_experiments

def setup_experiment(new_version_accessions: List[str], old_version_accessions: List[str]) -> Dict:
    """ Create an experiment where some samples were processed with the newest version of salmon and
    other with an older one.
    """
    # Create the experiment
    experiment_accession = 'SRP095529'
    data_dir = '/home/user/data_store/'
    experiment_dir = data_dir + experiment_accession
    experiment = Experiment.objects.create(
        accession_code=experiment_accession,
        technology='RNA-SEQ'
    )

    zebrafish = Organism.get_object_for_name("DANIO_RERIO")

    # Create the transcriptome processor and result:
    transcriptome_processor = Processor()
    transcriptome_processor.name = "Transcriptome"
    transcriptome_processor.version = "salmon 0.9.1"
    transcriptome_processor.docker_image = "dr_transcriptome"
    transcriptome_processor.environment = '{"some": "environment"}'
    transcriptome_processor.save()
    computational_result_short = ComputationalResult(processor=transcriptome_processor)
    computational_result_short.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = zebrafish
    organism_index.result = computational_result_short
    organism_index.absolute_directory_path = "/home/user/data_store/ZEBRAFISH_INDEX/SHORT"
    organism_index.salmon_version='salmon 0.9.1'
    organism_index.save()

    comp_file = ComputedFile()
    # This path will not be used because we already have the files extracted.
    comp_file.absolute_file_path = "/home/user/data_store/ZEBRAFISH_INDEX/SHORT/zebrafish_short.tar.gz"
    comp_file.result = computational_result_short
    comp_file.size_in_bytes=1337
    comp_file.sha1="ABC"
    comp_file.s3_key = "key"
    comp_file.s3_bucket = "bucket"
    comp_file.save()

    quant_processor = Processor()
    quant_processor.name = "Salmon Quant"
    quant_processor.version = "salmon 0.9.1"
    quant_processor.docker_image = "dr_salmon"
    quant_processor.environment = '{"some": "environment"}'
    quant_processor.save()

    for accession_code in old_version_accessions:
        sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database='SRA',
            technology='RNA-SEQ',
            platform_accession_code='IlluminaHiSeq1000'
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

        original_file = OriginalFile()
        original_file.filename=accession_code+'.SRA'
        original_file.source_filename=accession_code+'.SRA'
        original_file.save()

        OriginalFileSampleAssociation.objects.get_or_create(original_file=original_file, sample=sample)

        # Create and associate quant result and files.
        quant_result = ComputationalResult()
        quant_result.is_ccdl = True
        quant_result.processor = quant_processor
        quant_result.organism_index = organism_index # associate with OLD organism index
        quant_result.save()

        kv = ComputationalResultAnnotation()
        kv.data = {"index_length": "short"}
        kv.result = quant_result
        kv.is_public = True
        kv.save()

        # In prod the filename pattern will involve the timestamp
        # but here we're using the accession code so we can find
        # the archive file for the current sample.
        archive_filename = "result-" + accession_code + ".tar.gz"
        archive_file = ComputedFile()
        archive_file.filename = archive_filename
        archive_file.absolute_file_path = os.path.join(experiment_dir, archive_filename)
        archive_file.is_public = False
        archive_file.is_smashable = False
        archive_file.is_qc = False
        archive_file.result = quant_result
        archive_file.size_in_bytes = 12345
        archive_file.save()

        quant_file = ComputedFile()
        quant_file.filename = "quant.sf"
        quant_file.absolute_file_path = experiment_dir + "/quant_files/" + accession_code + "_output/quant.sf"
        quant_file.is_public = False
        quant_file.is_smashable = False
        quant_file.is_qc = False
        quant_file.result = quant_result
        quant_file.size_in_bytes = 12345
        quant_file.s3_bucket = "bucket"
        quant_file.s3_key = "key"
        quant_file.save()

        SampleResultAssociation.objects.get_or_create(
            sample=sample,
            result=quant_result
        )

    # Create another OrganismIndex with a newer version of
    transcriptome_processor = Processor()
    transcriptome_processor.name = "Transcriptome"
    transcriptome_processor.version = "salmon 0.13.1"
    transcriptome_processor.docker_image = "dr_transcriptome"
    transcriptome_processor.environment = '{"some": "environment"}'
    transcriptome_processor.save()
    computational_result_short = ComputationalResult(processor=transcriptome_processor)
    computational_result_short.save()

    organism_index = OrganismIndex()
    organism_index.index_type = "TRANSCRIPTOME_SHORT"
    organism_index.organism = zebrafish
    organism_index.result = computational_result_short
    organism_index.absolute_directory_path = "/home/user/data_store/ZEBRAFISH_INDEX/SHORT"
    organism_index.salmon_version='salmon 0.13.1' # DIFFERENT SALMON VERSION
    organism_index.save()

    comp_file = ComputedFile()
    # This path will not be used because we already have the files extracted.
    comp_file.absolute_file_path = "/home/user/data_store/ZEBRAFISH_INDEX/SHORT/zebrafish_short.tar.gz"
    comp_file.result = computational_result_short
    comp_file.size_in_bytes=1337
    comp_file.sha1="ABC"
    comp_file.s3_key = "key"
    comp_file.s3_bucket = "bucket"
    comp_file.save()

    for accession_code in new_version_accessions:
        sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database='SRA',
            technology='RNA-SEQ',
            platform_accession_code='IlluminaHiSeq1000'
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

        original_file = OriginalFile()
        original_file.filename=accession_code+'.SRA'
        original_file.source_filename=accession_code+'.SRA'
        original_file.save()

        OriginalFileSampleAssociation.objects.get_or_create(original_file=original_file, sample=sample)

        # Create and associate quant result and files.
        quant_result = ComputationalResult()
        quant_result.is_ccdl = True
        quant_result.processor = quant_processor
        quant_result.organism_index = organism_index # NEWER VERSION
        quant_result.save()

        kv = ComputationalResultAnnotation()
        kv.data = {"index_length": "short"}
        kv.result = quant_result
        kv.is_public = True
        kv.save()

        # In prod the filename pattern will involve the timestamp
        # but here we're using the accession code so we can find
        # the archive file for the current sample.
        archive_filename = "result-" + accession_code + ".tar.gz"
        archive_file = ComputedFile()
        archive_file.filename = archive_filename
        archive_file.absolute_file_path = os.path.join(experiment_dir, archive_filename)
        archive_file.is_public = False
        archive_file.is_smashable = False
        archive_file.is_qc = False
        archive_file.result = quant_result
        archive_file.size_in_bytes = 12345
        archive_file.save()

        quant_file = ComputedFile()
        quant_file.filename = "quant.sf"
        quant_file.absolute_file_path = experiment_dir + "/quant_files/" + accession_code + "_output/quant.sf"
        quant_file.is_public = False
        quant_file.is_smashable = False
        quant_file.is_qc = False
        quant_file.result = quant_result
        quant_file.size_in_bytes = 12345
        quant_file.s3_bucket = "bucket"
        quant_file.s3_key = "key"
        quant_file.save()

        SampleResultAssociation.objects.get_or_create(
            sample=sample,
            result=quant_result
        )

    return experiment

class RerunSalmonTestCase(TestCase):
    """
    Tests that new processor jobs are created for samples that belong to experiments that were
    processed with multiple versions of Salmon
    """
    def test_no_processor_job_needed(self):
        setup_experiment(['AA001', 'AA002'], [])
        update_salmon_all_experiments()

        # Verify that no jobs were created, because all samples had been processed with the latest version
        dl_jobs = DownloaderJob.objects.all()
        self.assertEqual(dl_jobs.count(), 0)

    def test(self):
        setup_experiment(['SS001'], ['SS002'])
        update_salmon_all_experiments()

        dl_jobs = DownloaderJob.objects.all()
        self.assertEqual(dl_jobs.count(), 1)

    def test_no_job_created_when_failed_job_exists(self):
        experiment = setup_experiment([], ['GSM001'])

        # create a failed job for that experiment
        processor_job = ProcessorJob()
        processor_job.pipeline_applied = ProcessorPipeline.SALMON
        processor_job.ram_amount = 1024
        processor_job.success = False
        processor_job.retried = False
        processor_job.no_retry = False
        processor_job.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = experiment.samples.first().original_files.first()
        assoc.processor_job = processor_job
        assoc.save()

        # Run command
        update_salmon_all_experiments()

        dl_jobs = DownloaderJob.objects.all()
        self.assertEqual(dl_jobs.count(), 0)

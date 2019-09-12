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
from data_refinery_foreman.foreman.management.commands import run_tximport

def run_tximport_at_progress_point(complete_accessions: List[str], incomplete_accessions: List[str]) -> Dict:
    """Create an experiment and associated objects and run tximport on it.

    Creates a sample for each accession contained in either input
    list. The samples in complete_accessions will be simlulated as
    already having salmon quant run on them. The samples in
    incomplete_accessions won't.
    """
    # Create the experiment
    experiment_accession = 'SRP095529'
    data_dir = '/home/user/data_store/'
    experiment_dir = data_dir + experiment_accession
    experiment = Experiment.objects.create(
        accession_code=experiment_accession,
        organism_names=['HOMO_SAPIENS'],
        technology='RNA-SEQ'
    )

    zebrafish = Organism.get_object_for_name("DANIO_RERIO")

    # Create the transcriptome processor and result:
    transcriptome_processor = Processor()
    transcriptome_processor.name = "Transcriptome"
    transcriptome_processor.version = "v9.9.9"
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
    organism_index.salmon_version = 'salmon 0.13.1'
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

    for accession_code in incomplete_accessions:
        sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database='SRA',
            technology='RNA-SEQ'
        )
        ExperimentSampleAssociation.objects.create(experiment=experiment, sample=sample)

        original_file = OriginalFile()
        original_file.filename=accession_code+'.SRA'
        original_file.source_filename=accession_code+'.SRA'
        original_file.save()

        OriginalFileSampleAssociation.objects.get_or_create(original_file=original_file, sample=sample)

    quant_processor = Processor()
    quant_processor.name = "Salmon Quant"
    quant_processor.version = "v9.9.9"
    quant_processor.docker_image = "dr_salmon"
    quant_processor.environment = '{"some": "environment"}'
    quant_processor.save()
    tximport_processor = Processor()
    tximport_processor.name = "Tximport"
    tximport_processor.version = "v9.9.9"
    tximport_processor.docker_image = "dr_salmon"
    tximport_processor.environment = '{"some": "environment"}'
    tximport_processor.save()

    # Create the already processed samples along with their
    # ComputationalResults and ComputedFiles. They don't need
    # original files for this test because we aren't going to run
    # salmon quant on them.
    for accession_code in complete_accessions:
        sample = Sample.objects.create(
            accession_code=accession_code,
            organism=zebrafish,
            source_database='SRA',
            technology='RNA-SEQ'
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

    # Setup is done, actually run the command.
    run_tximport.run_tximport()



class RunTximportTestCase(TestCase):
    """Tests that run_tximport only queues Tximport jobs for experiments which are ready.

    Some experiments are going to have samples that can't be
    processed. This means that tximport can't run on those unless we
    tell it to. We have a tximport only job to do this, but it should
    only be run on experiments with at least 25 samples where at least
    80% of the samples have been processed. We therefore run the
    tximport job on an experiment that is ready for it, one that has
    too few samples, and one that has too low of a copmpletion
    percent.
    """
    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.get_active_volumes')
    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.send_job')
    def test_early_tximport(self, mock_send_job, mock_get_active_volumes):
        """Tests that tximport jobs are created when the experiment is past the thresholds.

        Makes sure that when we should in fact run create tximport jobs that
        we do so, it works, and that it works even if there already is
        an existing tximport result for the experiment.

        So we run tximport on an experiment with 5 samples that aren't
        yet processed and 20 samples that are. By having 20/25 samples
        complete, we're just past both the numerical and percent
        cutoffs.
        """
        # First, set up our mocks to prevent network calls.
        mock_send_job.return_value = True
        active_volumes =  {"1", "2", "3"}
        mock_get_active_volumes.return_value = active_volumes

        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        pj = ProcessorJob.objects.all()[0]
        self.assertEqual(pj.pipeline_applied, ProcessorPipeline.TXIMPORT.value)

        # Verify that we attempted to send the jobs off to nomad
        mock_calls = mock_send_job.mock_calls
        self.assertEqual(len(mock_calls), 1)

        first_call_job_type = mock_calls[0][1][0]
        self.assertEqual(first_call_job_type, ProcessorPipeline.TXIMPORT)

    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.get_active_volumes')
    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.send_job')
    def test_tximport_percent_cutoff(self, mock_send_job, mock_get_active_volumes):
        """Tests logic for determining if tximport should be run early.

        This test is verifying that a tximport job won't be creatd if
        the experiment has too low of a percentage of processed
        samples to be elegible.

        So we run tximport on an experiment with 6 samples that aren't
        yet processed and 20 samples that are. By having 20/26 samples
        complete, we're past the numerical cutoff, but not past the
        percent cutoff.
        """
        # First, set up our mocks to prevent network calls.
        mock_send_job.return_value = True
        active_volumes =  {"1", "2", "3"}
        mock_get_active_volumes.return_value = active_volumes

        # Accessions SRR5125615-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125615",
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125621",
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertEqual(ProcessorJob.objects.all().count(), 0)

        # Verify that we didn't attempt to send the jobs off to nomad
        mock_calls = mock_send_job.mock_calls
        self.assertEqual(len(mock_calls), 0)

    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.get_active_volumes')
    @patch('data_refinery_foreman.foreman.management.commands.run_tximport.send_job')
    def test_tximport_numerical_cutoff(self, mock_send_job, mock_get_active_volumes):
        """Tests logic for determining if tximport should be run early.

        This test is verifying that a tximport job won't be created if
        the experiment has too few samples in it to be elegible.

        So we run tximport on an experiment with 5 samples that aren't
        yet processed and 19 that are. By having 19/24 samples
        complete, we're past the percent cutoff, but not past the
        numerical cutoff.
        """
        # First, set up our mocks to prevent network calls.
        mock_send_job.return_value = True
        active_volumes =  {"1", "2", "3"}
        mock_get_active_volumes.return_value = active_volumes

        # Accessions SRR5125616-SRR5125620 don't exist in SRA, but we
        # don't actually want to process them so it's okay.
        incomplete_accessions = [
            "SRR5125616",
            "SRR5125617",
            "SRR5125618",
            "SRR5125619",
            "SRR5125620",
        ]

        complete_accessions = [
            "SRR5125622",
            "SRR5125623",
            "SRR5125624",
            "SRR5125625",
            "SRR5125626",
            "SRR5125627",
            "SRR5125628",
            "SRR5125629",
            "SRR5125630",
            "SRR5125631",
            "SRR5125632",
            "SRR5125633",
            "SRR5125634",
            "SRR5125635",
            "SRR5125636",
            "SRR5125637",
            "SRR5125638",
            "SRR5125639",
            "SRR5125640",
        ]

        job_context = run_tximport_at_progress_point(complete_accessions, incomplete_accessions)

        # Confirm that this experiment is not ready for tximport yet,
        # because `salmon quant` is not run on 'fake_sample' and it
        # doens't have enough samples to have tximport run early.
        self.assertEqual(ProcessorJob.objects.all().count(), 0)

        # Verify that we didn't attempt to send the jobs off to nomad
        mock_calls = mock_send_job.mock_calls
        self.assertEqual(len(mock_calls), 0)

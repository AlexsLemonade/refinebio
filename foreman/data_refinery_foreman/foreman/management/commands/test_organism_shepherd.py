import json
import datetime

from django.core.management import call_command
from django.test import TransactionTestCase
from django.utils import timezone
from unittest.mock import MagicMock, Mock, patch, call

from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentOrganismAssociation,
    ExperimentResultAssociation,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    Organism,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
)

class OrganismShepherdTestCase(TransactionTestCase):
    @patch('data_refinery_foreman.foreman.management.commands.organism_shepherd.send_job')
    @patch('data_refinery_foreman.foreman.management.commands.organism_shepherd.Nomad')
    def test_organism_shepherd_command(self, mock_nomad, mock_send_job):
        """Tests that the organism shepherd requeues jobs in the right order.

        The situation we're setting up is basically this:
          * There are two experiments.
          * One of them has 1/2 samples processed, the other 0/1
          * One of them needs a DownloaderJob requeued and the other
            needs a ProcessorJob requued.

        And what we're going to test for is:
          * Both of the jobs that need to be requeued are requeued.
          * The experiment with a processed sample is requeued first
            because it has a higher completion percentage.
        """
        # First, set up our mocks to prevent network calls.
        mock_send_job.return_value = True

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.jobs = MagicMock()
            ret_value.jobs.get_jobs = MagicMock()
            ret_value.jobs.get_jobs.side_effect = lambda: []
            return ret_value

        mock_nomad.side_effect=mock_init_nomad
        zebrafish = Organism(name="DANIO_RERIO", taxonomy_id=1337, is_scientific_name=True)
        zebrafish.save()


        # Experiment that is 0% complete.
        zero_percent_experiment = Experiment(accession_code='ERP037000')
        zero_percent_experiment.technology = 'RNA-SEQ'
        zero_percent_experiment.save()

        organism_assoc = ExperimentOrganismAssociation.objects.create(
            organism=zebrafish,
            experiment=zero_percent_experiment
        )

        zero_percent = OriginalFile()
        zero_percent.filename = "ERR037001.fastq.gz"
        zero_percent.source_filename = "ERR037001.fastq.gz"
        zero_percent.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR037/ERR037001/ERR037001_1.fastq.gz"
        zero_percent.is_archive = True
        zero_percent.save()

        zero_percent_sample = Sample()
        zero_percent_sample.accession_code = 'ERR037001'
        zero_percent_sample.organism = zebrafish
        zero_percent_sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.sample = zero_percent_sample
        assoc.original_file = zero_percent
        assoc.save()

        assoc = ExperimentSampleAssociation()
        assoc.sample = zero_percent_sample
        assoc.experiment = zero_percent_experiment
        assoc.save()

        # TODO: fix names of all the variables to be appropriate for this test case.
        zero_percent_dl_job = DownloaderJob()
        zero_percent_dl_job.accession_code = zero_percent_sample.accession_code
        zero_percent_dl_job.downloader_task = "SRA"
        zero_percent_dl_job.start_time = timezone.now()
        zero_percent_dl_job.end_time = timezone.now()
        zero_percent_dl_job.success = False
        zero_percent_dl_job.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = zero_percent_dl_job
        assoc.original_file = zero_percent
        assoc.save()

        # Experiment that is 50% complete.
        fify_percent_experiment = Experiment(accession_code='ERP036000')
        fify_percent_experiment.technology = 'RNA-SEQ'
        fify_percent_experiment.save()

        organism_assoc = ExperimentOrganismAssociation.objects.create(
            organism=zebrafish,
            experiment=fify_percent_experiment
        )

        ## First sample, this one has been processed.
        successful_pj = ProcessorJob()
        successful_pj.accession_code = "ERR036000"
        successful_pj.pipeline_applied = "SALMON"
        successful_pj.ram_amount = 12288
        successful_pj.start_time = timezone.now()
        successful_pj.end_time = timezone.now()
        successful_pj.success = True
        successful_pj.save()

        successful_og = OriginalFile()
        successful_og.filename = "ERR036000.fastq.gz"
        successful_og.source_filename = "ERR036000.fastq.gz"
        successful_og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036000/ERR036000_1.fastq.gz"
        successful_og.is_archive = True
        successful_og.save()

        successful_sample = Sample()
        successful_sample.accession_code = 'ERR036000'
        successful_sample.organism = zebrafish
        successful_sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.sample = successful_sample
        assoc.original_file = successful_og
        assoc.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.processor_job = successful_pj
        assoc.original_file = successful_og
        assoc.save()

        assoc = ExperimentSampleAssociation()
        assoc.sample = successful_sample
        assoc.experiment = fify_percent_experiment
        assoc.save()

        ## Second sample, this one hasn't been processed.
        fifty_percent_unprocessed_og = OriginalFile()
        fifty_percent_unprocessed_og.filename = "ERR036001.fastq.gz"
        fifty_percent_unprocessed_og.source_filename = "ERR036001.fastq.gz"
        fifty_percent_unprocessed_og.source_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR036/ERR036001/ERR036001_1.fastq.gz"
        fifty_percent_unprocessed_og.is_archive = True
        fifty_percent_unprocessed_og.save()

        fifty_percent_unprocessed_sample = Sample()
        fifty_percent_unprocessed_sample.accession_code = 'ERR036001'
        fifty_percent_unprocessed_sample.organism = zebrafish
        fifty_percent_unprocessed_sample.save()

        assoc = OriginalFileSampleAssociation()
        assoc.sample = fifty_percent_unprocessed_sample
        assoc.original_file = fifty_percent_unprocessed_og
        assoc.save()

        assoc = ExperimentSampleAssociation()
        assoc.sample = fifty_percent_unprocessed_sample
        assoc.experiment = fify_percent_experiment
        assoc.save()

        fifty_percent_processor_job = ProcessorJob()
        fifty_percent_processor_job.pipeline_applied = "SALMON"
        fifty_percent_processor_job.accession_code = fifty_percent_unprocessed_sample.accession_code
        fifty_percent_processor_job.ram_amount = 12288
        fifty_percent_processor_job.start_time = timezone.now()
        fifty_percent_processor_job.end_time = timezone.now()
        fifty_percent_processor_job.success = False
        fifty_percent_processor_job.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.processor_job = fifty_percent_processor_job
        assoc.original_file = fifty_percent_unprocessed_og
        assoc.save()

        # Setup is done, actually run the command.
        args = []
        options = {"organism_name": "DANIO_RERIO"}
        call_command("organism_shepherd", *args, **options)

        # Verify that the jobs were called in the correct order.
        mock_calls = mock_send_job.mock_calls

        first_call_job_type = mock_calls[0][1][0]
        first_call_job_object = mock_calls[0][2]["job"]
        self.assertEqual(first_call_job_type, ProcessorPipeline.SALMON)
        self.assertEqual(first_call_job_object.pipeline_applied, fifty_percent_processor_job.pipeline_applied)
        self.assertEqual(first_call_job_object.ram_amount, fifty_percent_processor_job.ram_amount)

        fifty_percent_processor_job.refresh_from_db()
        self.assertEqual(first_call_job_object, fifty_percent_processor_job.retried_job)

        second_call_job_type = mock_calls[1][1][0]
        second_call_job_object = mock_calls[1][2]["job"]
        self.assertEqual(second_call_job_type, Downloaders.SRA)
        self.assertEqual(second_call_job_object.accession_code, zero_percent_dl_job.accession_code)
        self.assertEqual(second_call_job_object.downloader_task, zero_percent_dl_job.downloader_task)

        zero_percent_dl_job.refresh_from_db()
        self.assertEqual(second_call_job_object, zero_percent_dl_job.retried_job)

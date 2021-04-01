from time import sleep
from unittest import TestCase

from django.test import tag

from data_refinery_common.models import DownloaderJob, Experiment, ProcessorJob, Sample
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment

MICROARRAY_ACCESSION_CODES = [
    "E-TABM-496",  # 39 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE96849",  # 66 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE41094",  # 18 samples of SACCHAROMYCES_CEREVISIAE submitter processed data
]

RNA_SEQ_ACCESSION_CODES = [
    "SRP047410",  # 21 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
    "SRP094706",  # 4 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
]

EXPERIMENT_ACCESSION_CODES = MICROARRAY_ACCESSION_CODES + RNA_SEQ_ACCESSION_CODES


# Use unittest TestCase instead of django TestCase to avoid the test
# being done in a transaction.
class EndToEndTestCase(TestCase):
    """In order to parallelize the jobs as much as possible, everything is
    done in one big ol' function.

    This includes, in order:
      * Purging previously surveyed experiments.
      * Surveying experiments which will trigger:
        * An array_express downloader job triggering a affymetrix processor job.
        * A GEO downloader job triggering a affymetrix processor job.
        * A GEO downloader job triggering a NO_OP processor job.
      * Creating a transcriptome index for SACCHAROMYCES_CEREVISIAE.
      * Surveying experiments which will trigger:
        * A SRA downloader job triggering a salmon processor job for 21 samples.
          (TODO: fail one of the jobs and run a tximport job to finish the experiment.)
        * A SRA downloader job triggering a salmon processor job for 4 samples.
          (This should be fully processed without intervention.)
    TODO:
      * Creating a QN Target for SACCHAROMYCES_CEREVISIAE.
      * Creating a Compendium for SACCHAROMYCES_CEREVISIAE.
      * Creating a Quantpendium for SACCHAROMYCES_CEREVISIAE.
      * Downloading a dataset aggregated by species.
    """

    @tag("end_to_end")
    def test_all_the_things(self):
        for accession_code in EXPERIMENT_ACCESSION_CODES:
            purge_experiment(accession_code)

        survey_jobs = []
        # Kick off microarray jobs first because downloading the
        # affymetrix image takes a long time.
        for accession_code in MICROARRAY_ACCESSION_CODES:
            survey_jobs.append(queue_surveyor_for_accession(accession_code))

        # However, before we kick of RNA-Seq jobs, we have to build a transcriptome index for them.
        transcriptome_survey_job = queue_surveyor_for_accession(
            "SACCHAROMYCES_CEREVISIAE, EnsemblFungi"
        )

        while transcriptome_survey_job.retried_job is None:
            sleep(20)
            print("Polling original transcriptome survey job.")
            transcriptome_survey_job.refresh_from_db()

        retried_transcritpome_survey_job = transcriptome_survey_job.retried_job

        while retried_transcritpome_survey_job.end_time is None:
            sleep(20)
            print("Polling retried transcriptome survey job.")
            retried_transcritpome_survey_job.refresh_from_db()

        transcriptome_downloader_jobs = DownloaderJob.objects.filter(
            downloader_task="TRANSCRIPTOME_INDEX",
            created_at__gt=retried_transcritpome_survey_job.created_at,
        )

        for downloader_job in transcriptome_downloader_jobs:
            while downloader_job.end_time is None:
                sleep(20)
                print("Polling transcriptome downloader jobs.")
                downloader_job.refresh_from_db()

        transcriptome_processor_jobs = ProcessorJob.objects.filter(
            processor_pipeline__startswith="TRANSCRIPTOME_INDEX",
            created_at__gt=retried_transcritpome_survey_job.created_at,
        )

        for processor_job in transcriptome_processor_jobs:
            while processor_job.end_time is None:
                sleep(20)
                print("Polling processor jobs.")
                processor_job.refresh_from_db()

        survey_jobs = []
        for accession_code in RNA_SEQ_ACCESSION_CODES:
            survey_jobs.append(queue_surveyor_for_accession(accession_code))

        retried_jobs = []
        for survey_job in survey_jobs:
            while survey_job.retried_job is None:
                sleep(20)
                print("Polling original survey jobs.")
                survey_job.refresh_from_db()

            retried_jobs.append(survey_job.retried_job)

        for retried_job in retried_jobs:
            while retried_job.end_time is None:
                sleep(20)
                print("Polling retried survey jobs.")
                retried_job.refresh_from_db()

        self.assertEqual(Sample.objects.count(), 147)

        samples = []
        for accession_code in EXPERIMENT_ACCESSION_CODES:
            experiment = Experiment.objects.get(accession_code=accession_code)
            samples.extend(list(experiment.samples.all()))

        downloader_jobs = []
        for sample in samples:
            greatest_job_id = -1
            last_job = None

            # There should be only one, but if it got retried for some
            # reason we want the latest one.
            for job in sample.get_downloader_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            downloader_jobs.append(last_job)

        for downloader_job in downloader_jobs:
            while downloader_job.end_time is None:
                sleep(20)
                print("Polling downloader jobs.")
                downloader_job.refresh_from_db()

        processor_jobs = []
        for sample in samples:
            greatest_job_id = -1
            last_job = None

            for job in sample.get_processor_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            processor_jobs.append(last_job)

        for processor_job in processor_jobs:
            while processor_job.end_time is None:
                sleep(20)
                print("Polling processor jobs.")
                processor_job.refresh_from_db()

        self.assertTrue(Sample.processed_objects.count() == 147)

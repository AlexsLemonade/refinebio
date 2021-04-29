from time import sleep
from unittest import TestCase

from django.test import tag

from data_refinery_common.models import DownloaderJob, Experiment, ProcessorJob, Sample
from data_refinery_foreman.foreman.management.commands.run_tximport import run_tximport
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment

MICROARRAY_ACCESSION_CODES = [
    "E-TABM-496",  # 39 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE96849",  # 68 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE41094",  # 18 samples of SACCHAROMYCES_CEREVISIAE submitter processed data
]

RNA_SEQ_ACCESSION_CODES = [
    "SRP047410",  # 26 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data, one will fail.
    "SRP094706",  # 4 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
]

EXPERIMENT_ACCESSION_CODES = MICROARRAY_ACCESSION_CODES + RNA_SEQ_ACCESSION_CODES


def wait_for_job(job) -> bool:
    """Waits for a job and all of its retries."""
    job.refresh_from_db()
    is_done = False

    while not is_done:

        if job.end_time is None and job.success is None:
            print(f"Polling {type(job).__name__}s. Currently waiting for job id: {job.id}")
            sleep(20)
            job.refresh_from_db()
        elif job.retried and job.retried_job:
            job = job.retried_job
        elif job.success:
            return True
        else:
            return False

    return False


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
        transcriptome_survey_job = queue_surveyor_for_accession("SACCHAROMYCES_CEREVISIAE, Ensembl")

        print(
            "First, creating transcriptome indices (and starting on Affy jobs while we're at it):"
        )
        self.assertTrue(wait_for_job(transcriptome_survey_job))

        transcriptome_downloader_jobs = DownloaderJob.objects.filter(
            downloader_task="TRANSCRIPTOME_INDEX",
            created_at__gt=transcriptome_survey_job.created_at,
        )

        for downloader_job in transcriptome_downloader_jobs:
            self.assertTrue(wait_for_job(downloader_job))

        transcriptome_processor_jobs = ProcessorJob.objects.filter(
            pipeline_applied__startswith="TRANSCRIPTOME_INDEX",
            created_at__gt=transcriptome_survey_job.created_at,
        )

        for processor_job in transcriptome_processor_jobs:
            self.assertTrue(wait_for_job(processor_job))

        for accession_code in RNA_SEQ_ACCESSION_CODES:
            survey_jobs.append(queue_surveyor_for_accession(accession_code))

        print("Next, processing all the raw data.")
        for survey_job in survey_jobs:
            self.assertTrue(wait_for_job(survey_job))

        self.assertEqual(Sample.objects.count(), 155)

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
            self.assertTrue(wait_for_job(downloader_job))

        processor_jobs = []
        for sample in samples:
            # This sample fails for good reason, so don't expect it to pass.
            if sample.accession_code == "SRR1583739":
                continue

            greatest_job_id = -1
            last_job = None

            for job in sample.get_processor_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            processor_jobs.append(last_job)

        for processor_job in processor_jobs:
            self.assertTrue(wait_for_job(processor_job))

        # Because SRR1583739 fails, the 26 samples from SRP047410 won't be processed
        self.assertEqual(Sample.processed_objects.count(), 129)

        print("Finally, need to run tximport to finish an experiment with one bad sample.")
        tximport_jobs = run_tximport()
        self.assertEqual(len(tximport_jobs), 1)

        self.assertTrue(wait_for_job(tximport_jobs[0]))

        # This is the full total of jobs minus one.
        self.assertEqual(Sample.processed_objects.count(), 154)

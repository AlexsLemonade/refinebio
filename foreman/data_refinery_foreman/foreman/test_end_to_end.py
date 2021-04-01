from time import sleep
from unittest import TestCase

from django.test import tag

from data_refinery_common.models import Sample
from data_refinery_foreman.surveyor.management.commands.surveyor_dispatcher import (
    queue_surveyor_for_accession,
)
from data_refinery_foreman.surveyor.management.commands.unsurvey import purge_experiment

EXPERIMENT_ACCESSIONS = [
    "E-TABM-496",  # 39 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "GSE96849",  # 66 samples of SACCHAROMYCES_CEREVISIAE microarray data
    "SRP076395",  # 20 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
    "SRP094706",  # 4 samples of SACCHAROMYCES_CEREVISIAE RNA-Seq data
    "GSE41094",  # 18 samples of SACCHAROMYCES_CEREVISIAE submitter processed data
]


# Use unittest TestCase instead of django TestCase to avoid the test
# being done in a transaction.
class EndToEndTestCase(TestCase):
    @tag("end_to_end")
    def test_testing(self):
        purge_experiment("GSE69775")
        survey_job = queue_surveyor_for_accession("GSE69775")

        while survey_job.retried_job is None:
            sleep(20)
            print("Polling original survey job.")
            survey_job.refresh_from_db()

        retried_job = survey_job.retried_job

        while retried_job.end_time is None:
            sleep(20)
            print("Polling retried survey job.")
            retried_job.refresh_from_db()

        self.assertTrue(Sample.objects.count() == 2)

        downloader_jobs = []
        for sample in Sample.objects.all():
            greatest_job_id = -1
            last_job = None

            for job in sample.get_downloader_jobs():
                if job.id > greatest_job_id:
                    greatest_job_id = job.id
                    last_job = job

            downloader_jobs.append(last_job)

        for downloader_job in downloader_jobs:
            while downloader_job.end_time is None:
                sleep(20)
                print("Polling downloader job.")
                downloader_job.refresh_from_db()

        processor_jobs = []
        for sample in Sample.objects.all():
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
                print("Polling processor job.")
                processor_job.refresh_from_db()

        self.assertTrue(Sample.processed_objects.count() == 2)

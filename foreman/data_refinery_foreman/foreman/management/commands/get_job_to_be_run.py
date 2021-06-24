import json

from django.core.management.base import BaseCommand

from data_refinery_common.models import DownloaderJob, ProcessorJob


class Command(BaseCommand):
    def handle(self, *args, **options):
        """Outputs a JSON dict of a job that need to be run.

        The most recently created unstarted ProcessorJob will be prioritized.

        If there are no unstarted ProcessorJobs, then the most
        recently created unstarted DownloaderJob will be prioritized.

        The JSON output will be a dict containing the following keys:
        * job_name
        * job_id
        * image_name
        * job_type
        """
        processor_jobs_to_run = [
            "AFFY_TO_PCL",
            "SALMON",
            "TXIMPORT",
            "ILLUMINA_TO_PCL",
            "TRANSCRIPTOME_INDEX_LONG",
            "TRANSCRIPTOME_INDEX_SHORT",
            "NO_OP",
        ]
        processor_job = (
            ProcessorJob.objects.filter(
                start_time=None, success=None, pipeline_applied__in=processor_jobs_to_run
            )
            .order_by("-id")
            .first()
        )
        if processor_job:
            job_dict = {
                "job_id": processor_job.id,
                "job_type": "ProcessorJob",
                "job_name": processor_job.pipeline_applied,
            }
            print(json.dumps(job_dict))
            return

        downloader_job = (
            DownloaderJob.objects.filter(start_time=None, success=None).order_by("-id").first()
        )
        if downloader_job:
            job_dict = {
                "job_id": downloader_job.id,
                "job_type": "DownloaderJob",
                "job_name": downloader_job.downloader_task,
            }
            print(json.dumps(job_dict))
            return

        print("{}")

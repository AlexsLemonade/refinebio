from django.core.management.base import BaseCommand

from data_refinery_common.message_queue import terminate_job
from data_refinery_common.models import DownloaderJob, ProcessorJob, SurveyJob

PROCESSOR = "PROCESSOR"
DOWNLOADER = "DOWNLOADER"
SURVEYOR = "SURVEYOR"


def stop_job(job_id, job_class, reason):
    job = job_class.objects.get(pk=job_id)

    if job.success is not None:
        # The job already completed, this would just further confuse things.
        return

    job.no_retry = True
    job.success = False
    job.failure_reason = reason
    job.save()

    terminate_job(job, reason)


def stop_jobs(job_ids, job_type, reason):
    if job_type == PROCESSOR:
        job_class = ProcessorJob
    elif job_type == DOWNLOADER:
        job_class = DownloaderJob
    elif job_type == SURVEYOR:
        job_class = SurveyJob
    else:
        print("Invalid job-type supplied!")
        exit(1)

    for job_id in job_ids.strip().split(","):
        stop_job(job_id, job_class, reason)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "job-ids", type=str, help=("Comma separated job IDs that need to be stopped."),
        )

        parser.add_argument(
            "--job-type",
            type=str,
            default=PROCESSOR,
            choices=[PROCESSOR, DOWNLOADER, SURVEYOR],
            help=(
                f"The type of job to be stopped, either {PROCESSOR}, {DOWNLOADER}, or {SURVEYOR}."
            ),
        )

        parser.add_argument(
            "--reason",
            type=str,
            default=(
                "The job was terminated by"
                " data_refinery_foreman.foreman.management.commands.stop_job"
            ),
            help=(
                "A reason to be used as the job's failure_reason and to be provided to"
                " AWS Batch as the reason for the job termination."
            ),
        )

    def handle(self, *args, **options):
        """This command stops the job from running and marks it as no retry.

        Defaults to stopping processor jobs,

        If a comma-separated --accesion-codes parameter is supplied
        only eligible experiments in that list will be run instead.
        """
        stop_jobs(options["job-ids"], options["job_type"], options["reason"])

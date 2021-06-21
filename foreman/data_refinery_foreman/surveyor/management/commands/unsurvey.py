"""
This command will remove experiments and all data that is uniquely associated with them.
(i.e. it's only associated with these experiment and not others.)
This can be used so that a surveyor job that went wrong can be rerun.
"""
import sys
import uuid

from django.core.management.base import BaseCommand

import boto3
import botocore

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import *
from data_refinery_common.utils import parse_s3_url

logger = get_and_configure_logger(__name__)


def delete_job_and_retries(job) -> None:
    """Deletes a job and any jobs that retried it."""
    try:
        retried_job = job.retried_job
    except Exception:
        # I'm not sure why this isn't safe. It might be that we're
        # somehow able to delete the job this is referencing.
        retried_job = None

    try:
        job.delete()
    except Exception:
        # If the job already was deleted it's quite okay.
        pass

    if retried_job:
        delete_job_and_retries(retried_job)


def purge_experiment(accession: str) -> None:
    """Removes an experiment and all data that is only associated with it.

    I.e. if a sample is part of two experiments, it won't be removed
    but it will be disassociated from the one being purged.
    """
    experiment = Experiment.objects.filter(accession_code=accession)[0]
    ExperimentAnnotation.objects.filter(experiment=experiment).delete()

    ## Samples
    experiment_sample_assocs = ExperimentSampleAssociation.objects.filter(experiment=experiment)
    samples = Sample.objects.filter(id__in=experiment_sample_assocs.values("sample_id"))

    uniquely_assoced_sample_ids = []
    for sample in samples:
        assoc_query = ExperimentSampleAssociation.objects.filter(sample=sample)
        if assoc_query.count() == 1:
            # The sample is only associated with this experiment, mark it for deletion.
            uniquely_assoced_sample_ids.append(sample.id)
            SampleAnnotation.objects.filter(sample=sample).delete()
            assoc_query[0].delete()
        else:
            # The association isn't unique, but we still should delete
            # the association with the experiment since it will be deleted.
            ExperimentSampleAssociation.objects.filter(sample=sample, experiment=experiment)[
                0
            ].delete()

    ## ComputationalResults/ComputedFiles
    comp_result_sample_assocs = SampleResultAssociation.objects.filter(
        sample_id__in=uniquely_assoced_sample_ids
    )
    comp_results = ComputationalResult.objects.filter(
        id__in=comp_result_sample_assocs.values("result_id")
    )

    uniquely_assoced_comp_result_ids = []
    for comp_result in comp_results:
        # We know this result is associated with samples we can
        # delete, but is it associated with any we cannot?
        extra_assocs = SampleResultAssociation.objects.filter(result=comp_result).exclude(
            sample_id__in=uniquely_assoced_sample_ids
        )
        if extra_assocs.count() == 0:
            # It's not associated with anything else, delete it, its
            # associations, and its ComputedFile.
            uniquely_assoced_comp_result_ids.append(comp_result.id)

            ComputationalResultAnnotation.objects.filter(result=comp_result).delete()
            computed_files = ComputedFile.objects.filter(result=comp_result)

            # Delete all of the local files
            for computed_file in computed_files:
                computed_file.delete_local_file(force=True)

            # Delete the database records for the ComputedFile
            SampleComputedFileAssociation.objects.filter(computed_file__in=computed_files).delete()

            # pre_delete will also delete the files on s3
            computed_files.delete()

    # Whether or not we can delete all of these results, we know the
    # associations need to go.
    comp_result_sample_assocs.delete()

    # OriginalFiles
    og_file_sample_assocs = OriginalFileSampleAssociation.objects.filter(
        sample_id__in=uniquely_assoced_sample_ids
    )
    original_files = OriginalFile.objects.filter(
        id__in=og_file_sample_assocs.values("original_file_id")
    )

    uniquely_assoced_og_file_ids = []
    for original_file in original_files:
        # We know this file is associated with samples we can
        # delete, but is it associated with any we cannot?
        extra_assocs = OriginalFileSampleAssociation.objects.filter(
            original_file=original_file
        ).exclude(sample_id__in=uniquely_assoced_sample_ids)
        if extra_assocs.count() == 0:
            # Build a list of original_files so we can delete them all at once.
            uniquely_assoced_og_file_ids.append(original_file.id)

            # However we have to delete their local files one at a time:
            try:
                # Not all original files have absolute_file_path set apparently
                original_file.delete_local_file()
            except Exception:
                pass

    # Whether or not we can delete all of these original files, we
    # know the associations need to go.
    og_file_sample_assocs.delete()

    ## DownloaderJobs
    og_file_dj_assocs = DownloaderJobOriginalFileAssociation.objects.filter(
        original_file_id__in=uniquely_assoced_og_file_ids
    )
    # Important to order by id, so the jobs that didn't retry anything are deleted first.
    downloader_jobs = DownloaderJob.objects.filter(
        id__in=og_file_dj_assocs.values("downloader_job_id")
    ).order_by("id")

    for downloader_job in downloader_jobs:
        # We know this job is associated with files we can
        # delete, but is it associated with any we cannot?
        extra_assocs = DownloaderJobOriginalFileAssociation.objects.filter(
            downloader_job=downloader_job
        ).exclude(original_file_id__in=uniquely_assoced_og_file_ids)
        if extra_assocs.count() == 0:
            delete_job_and_retries(downloader_job)

    # Whether or not we can delete all of these jobs, we
    # know the associations need to go.
    og_file_dj_assocs.delete()

    ## ProcessorJobs
    og_file_pj_assocs = ProcessorJobOriginalFileAssociation.objects.filter(
        original_file_id__in=uniquely_assoced_og_file_ids
    )
    # Important to order by id, so the jobs that didn't retry anything are deleted first.
    processor_jobs = ProcessorJob.objects.filter(
        id__in=og_file_pj_assocs.values("processor_job_id")
    ).order_by("id")

    for processor_job in processor_jobs:
        # We know this job is associated with files we can
        # delete, but is it associated with any we cannot?
        extra_assocs = ProcessorJobOriginalFileAssociation.objects.filter(
            processor_job=processor_job
        ).exclude(original_file_id__in=uniquely_assoced_og_file_ids)
        if extra_assocs.count() == 0:
            delete_job_and_retries(processor_job)

    # Whether or not we can delete all of these jobs, we
    # know the associations need to go.
    og_file_pj_assocs.delete()

    # Delete the lists of objects we built.
    OriginalFile.objects.filter(id__in=uniquely_assoced_og_file_ids).delete()
    ComputationalResult.objects.filter(id__in=uniquely_assoced_comp_result_ids).delete()
    Sample.objects.filter(id__in=uniquely_assoced_sample_ids).delete()

    # Finally delete the experiment itself last.
    experiment.delete()

    logger.info("Purged experiments: %s", accession)


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession", help=("An experiment accession code to survey, download, and process.")
        )
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes. s3:// URLs are also accepted."),
        )
        parser.add_argument("--force", help="Do not ask for confirmation.", action="store_true")

    def handle(self, *args, **options):
        if options["accession"] is None and options["file"] is None:
            logger.error("You must specify an experiment accession or file.")
            sys.exit(1)

        if not options["force"]:
            print("-------------------------------------------------------------------------------")
            print(
                "This will delete all objects in the database related to these accessions."
                " Are you sure you want to do this?"
            )
            answer = input('You must type "yes", all other input will be ignored: ')

            if answer != "yes":
                print("Not unsurveying because confirmation was denied.")
                sys.exit(1)

        accessions = []
        if options["file"]:
            if "s3://" in options["file"]:
                bucket, key = parse_s3_url(options["file"])
                s3 = boto3.resource("s3")
                try:
                    filepath = "/tmp/input_" + str(uuid.uuid4()) + ".txt"
                    s3.Bucket(bucket).download_file(key, filepath)
                except botocore.exceptions.ClientError as e:
                    if e.response["Error"]["Code"] == "404":
                        logger.error("The remote file does not exist.")
                        raise
                    else:
                        raise
            else:
                filepath = options["file"]

            with open(filepath) as file:
                for accession in file:
                    accessions.append(accession.strip())
        else:
            accessions.append(options["accession"])

        for accession in accessions:
            logger.info("Purging Experiment with accession: %s", accession)
            try:
                purge_experiment(accession)
            except Exception:
                logger.exception(
                    "Exception caught while purging experiment with accession: %s", accession
                )

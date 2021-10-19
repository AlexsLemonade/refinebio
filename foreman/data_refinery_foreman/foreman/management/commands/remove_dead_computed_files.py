import time

from django.core.management.base import BaseCommand
from django.db.models import Q

from data_refinery_common.job_management import create_downloader_job
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.performant_pagination.pagination import PAGE_SIZE, PerformantPaginator


def requeue_sample(sample, dry_run=False):
    sample_was_requeued = False
    if sample.is_processed:
        has_live_computational_results = False
        for result in sample.results.all():
            live_files = result.computedfile_set.filter(
                s3_bucket__isnull=False, s3_key__isnull=False
            )
            if live_files.count() >= 1:
                has_live_computational_results = True

        live_computed_files = sample.computed_files.filter(
            s3_bucket__isnull=False, s3_key__isnull=False
        )

        if not (has_live_computational_results or live_computed_files.count() > 1):
            sample_was_requeued = True
            if not dry_run:
                # There's no live computed files, the sample
                # should not have been marked processed.
                sample.is_processed = False
                sample.save()

                create_downloader_job(sample.original_files.all(), force=True)

    return sample_was_requeued


def requeue_samples(sample_queryset, dry_run=False):
    paginator = PerformantPaginator(sample_queryset, PAGE_SIZE)
    page = paginator.page()

    # Loop through the samples one by one to see if they've been
    # erroneously marked as processed. If so, mark them as
    # unprocessed, and kick off a new job so they can get
    # processed correctly.
    # Do this before deleting the computed files in case we get
    # interrupted. It'll be harder to tell what samples were
    # erroneously marked as processed.
    while True:
        counter = 0
        for sample in page.object_list:
            if requeue_sample(sample, dry_run):
                counter += 1
            # requeue_sample makes database calls, not a good idea to
            # call in a loop without a sleep.
            time.sleep(1)

        print(f"Requeued {counter} samples in that page.")

        if not page.has_next():
            break
        else:
            page = paginator.page(page.next_page_number())


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            help="Prints what it would do, without doing any of it.",
            action="store_true",
        )

    def handle(self, *args, **options):
        """Removes computed files that weren't uploaded to S3.

        It also cleans up computational results that are relying on these computed files.

        It also marks samples that are processed as unprocessed if
        they don't have another live computed file or result, and then
        requeues a downloader job for them to be processed and have
        their files uploaded to s3.
        """
        dead_computed_files = ComputedFile.objects.filter(
            Q(s3_bucket__isnull=True) | Q(s3_key__isnull=True)
        )

        file_sample_assocs = SampleComputedFileAssociation.objects.filter(
            computed_file_id__in=dead_computed_files.values("id")
        )
        file_associated_samples = Sample.objects.filter(
            id__in=file_sample_assocs.values("sample_id")
        )

        requeue_samples(file_associated_samples, options["dry_run"])

        # TODO: Clear out ComputationalResults and their associated compendium results.
        # What about file->result->sample? Don't some of thems not get associated directly?

        # Some ComputedFiles are only linked to samples indirectly
        # through their ComputationalResults. We need also need to run
        # requeue_sample on them.
        dead_computational_results = ComputationalResult.objects.filter(
            id__in=dead_computed_files.values("result_id").distinct()
        )
        result_sample_assocs = SampleResultAssociation.objects.filter(
            result_id__in=dead_computational_results.values("id")
        )
        result_associated_samples = Sample.objects.filter(
            id__in=result_sample_assocs.values("sample_id")
        )

        requeue_samples(result_associated_samples, options["dry_run"])

        print(f"Deleting {dead_computational_results.count()} ComputationalResults")
        print(f"Deleting {dead_computed_files.count()} ComputedFiles")
        if not options["dry_run"]:
            dead_computational_results.delete()
            dead_computed_files.delete()

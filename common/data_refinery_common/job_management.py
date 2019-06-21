from typing import List, Dict, Callable

from data_refinery_common import job_lookup
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
)
from data_refinery_common.utils import (
    get_env_variable,
    get_env_variable_gracefully,
    get_instance_id,
)


logger = get_and_configure_logger(__name__)


def create_downloader_job(undownloaded_files: List[OriginalFile],
                          *,
                          processor_job_id=None,
                          force=False) -> bool:
    """Creates a downloader job to download `undownloaded_files`."""
    if not undownloaded_files:
        return False

    original_downloader_job = None
    archive_file = None
    for undownloaded_file in undownloaded_files:
        try:
            original_downloader_job = undownloaded_file.downloader_jobs.latest('id')

            # Found the job so we don't need to keep going.
            break
        except DownloaderJob.DoesNotExist:
            # If there's no association between this file and any
            # downloader jobs, it's most likely because the original
            # file was created after extracting a archive containing
            # multiple files worth of data.
            # The way to handle this is to find that archive and
            # recreate a downloader job FOR THAT. That archive will
            # have the same filename as the file at the end of the
            # 'source_url' field, because that source URL is pointing
            # to the archive we need.
            archive_filename = undownloaded_file.source_url.split("/")[-1]

            # This file or its job might not exist, but we'll wait
            # until we've checked all the files before calling it a
            # failure.
            try:
                archive_file = OriginalFile.objects.filter(filename=archive_filename)
                if archive_file.count() > 0:
                    archive_file = archive_file.first()
                else:
                    # We might need to match these up based on
                    # source_filenames rather than filenames so just
                    # try them both.
                    archive_file = OriginalFile.objects.filter(source_filename=archive_filename).first()

                original_downloader_job = DownloaderJobOriginalFileAssociation.objects.filter(
                    original_file=archive_file
                ).latest('id').downloader_job
                # Found the job so we don't need to keep going.
                break
            except:
                pass

    if not original_downloader_job:
        sample_object = list(undownloaded_files)[0].samples.first()
        if sample_object:
            downloader_task = job_lookup.determine_downloader_task(sample_object)

            if downloader_task == job_lookup.Downloaders.NONE:
                logger.warn(("No valid downloader task found for sample, which is weird"
                             " because it was able to have a processor job created for it..."),
                            sample=sample_object.id)
                return False
            else:
                # determine_downloader_task returns an enum object,
                # but we wanna set this on the DownloaderJob object so
                # we want the actual value.
                downloader_task = downloader_task.value

            accession_code = sample_object.accession_code
            original_files = sample_object.original_files.all()
        else:
            logger.error(
                "Could not find the original DownloaderJob or Sample for these files.",
                undownloaded_file=undownloaded_files
            )
            return False
    elif original_downloader_job.was_recreated and not force:
        logger.warn(
            "Downloader job has already been recreated once, not doing it again.",
            original_downloader_job=original_downloader_job,
            undownloaded_files=undownloaded_files
        )
        return False
    else:
        downloader_task = original_downloader_job.downloader_task
        accession_code = original_downloader_job.accession_code
        original_files = original_downloader_job.original_files.all()

        sample_object = original_files[0].samples.first()

    new_job = DownloaderJob()
    new_job.downloader_task = downloader_task
    new_job.accession_code = accession_code
    new_job.was_recreated = True
    new_job.ram_amount = 1024
    new_job.save()

    if archive_file:
        # If this downloader job is for an archive file, then the
        # files that were passed into this function aren't what need
        # to be directly downloaded, they were extracted out of this
        # archive. The DownloaderJob will re-extract them and set up
        # the associations for the new ProcessorJob.
        # So double check that it still needs downloading because
        # another file that came out of it could have already
        # recreated the DownloaderJob.
        if archive_file.needs_downloading(processor_job_id):
            if archive_file.is_downloaded:
                # If it needs to be downloaded then it's not
                # downloaded and the is_downloaded field should stop
                # lying about that.
                archive_file.is_downloaded = False
                archive_file.save()

            DownloaderJobOriginalFileAssociation.objects.get_or_create(
                downloader_job=new_job,
                original_file=archive_file
            )
    else:
        # We can't just associate the undownloaded files, because
        # there's a chance that there is a file which actually is
        # downloaded that also needs to be associated with the job.
        for original_file in original_files:
            DownloaderJobOriginalFileAssociation.objects.get_or_create(
                downloader_job=new_job,
                original_file=original_file
            )

    return True

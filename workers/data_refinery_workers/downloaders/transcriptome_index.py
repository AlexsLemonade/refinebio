from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
from typing import List
from contextlib import closing
from data_refinery_common.models import File, DownloaderJob
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _verify_files(file1: File, file2: File, job: DownloaderJob) -> None:
    """Verifies that the two files are the same.

    This is useful for this downloader because each job has two
    batches which should each have the same two files.
    """
    if file1.download_url != file2.download_url:
        failure_message = ("A Batch's file doesn't have the same download "
                           "URL as the other batch's file.")
        logger.error(failure_message,
                     downloader_job=job.id)
        job.failure_reason = failure_message
        raise ValueError(failure_message)


def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    failure_template = "Exception caught while downloading file from: %s"
    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path,
                     downloader_job=job.id)
        urllib.request.urlcleanup()
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()
    except Exception:
        logger.exception(failure_template,
                         download_url,
                         downloader_job=job.id)
        job.failure_reason = failure_template % download_url
        raise
    finally:
        target_file.close()


def _upload_files(job_dir: str, files: List[File], job: DownloaderJob) -> None:
    try:
        for file in files:
            file.size_in_bytes = os.path.getsize(file.get_temp_pre_path(job_dir))
            file.save()
            file.upload_raw_file(job_dir)
    except Exception:
        logger.exception("Exception caught while uploading file.",
                         downloader_job=job.id,
                         batch=file.batch.id)
        job.failure_reason = "Exception caught while uploading file."
        raise
    finally:
        file.remove_temp_directory(job_dir)


def download_transcriptome(job_id: int) -> None:
    """The main function for the Transcriptome Index Downloader.

    Two files are needed for the Transcriptome Index Downloader: a
    fasta file and a gtf file. However each pair need to be processed
    into two different sized indices. (See the
    processors.transcriptome_index._create_index function's docstring
    for more info.) Therefore we only download each set once, then
    push it to Temporary Storage twice.
    """
    job = utils.start_job(job_id)
    batches = job.batches.all()
    success = True
    job_dir = utils.JOB_DIR_PREFIX + str(job_id)

    try:
        first_fasta_file = File.objects.get(batch=batches[0], raw_format__exact="fa.gz")
        first_gtf_file = File.objects.get(batch=batches[0], raw_format__exact="gtf.gz")
        second_fasta_file = File.objects.get(batch=batches[1], raw_format__exact="fa.gz")
        second_gtf_file = File.objects.get(batch=batches[1], raw_format__exact="gtf.gz")
        os.makedirs(first_fasta_file.get_temp_dir(job_dir), exist_ok=True)
    except Exception:
        logger.exception("Failed to retrieve all expected files from database.",
                         downloader_job=job.id)
        job.failure_reason = "Failed to retrieve all expected files from database."
        success = False

    if success:
        try:
            _verify_files(first_fasta_file, second_fasta_file, job)
            _verify_files(first_gtf_file, second_gtf_file, job)

            # The two Batches share the same fasta and gtf files, so
            # only download each one once
            _download_file(first_fasta_file.download_url,
                           first_fasta_file.get_temp_pre_path(job_dir),
                           job)
            _download_file(first_gtf_file.download_url,
                           first_gtf_file.get_temp_pre_path(job_dir),
                           job)

            # Then create symlinks so the files for the second Batch
            # can be found where they will be expected to.
            try:
                os.symlink(first_fasta_file.get_temp_pre_path(job_dir),
                           second_fasta_file.get_temp_pre_path(job_dir))
                os.symlink(first_gtf_file.get_temp_pre_path(job_dir),
                           second_gtf_file.get_temp_pre_path(job_dir))
            except Exception:
                logger.exception("Exception caught while creating symlinks.",
                                 downloader_job=job.id)
                job.failure_reason = "Exception caught while creating symlinks."
                raise

            _upload_files(job_dir,
                          [first_fasta_file, first_gtf_file, second_fasta_file, second_gtf_file],
                          job)
        except Exception:
            # Exceptions are already logged and handled.
            # Just need to mark the job as failed.
            success = False

    if success:
        logger.debug("Files %s and %s downloaded successfully.",
                     first_fasta_file,
                     first_gtf_file,
                     downloader_job=job_id)

    utils.end_job(job, batches, success)

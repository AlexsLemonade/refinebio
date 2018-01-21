from __future__ import absolute_import, unicode_literals
import urllib.request
import os
import shutil
import zipfile
from typing import List
from contextlib import closing
from data_refinery_common.models import File, DownloaderJob
from data_refinery_workers.downloaders import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _verify_batch_grouping(files: List[File], job: DownloaderJob) -> None:
    """All batches in the same job should have the same downloader url"""
    for file in files:
        if file.download_url != files[0].download_url:
            failure_message = ("A Batch's file doesn't have the same download "
                               "URL as the other batches' files.")
            logger.error(failure_message,
                         downloader_job=job.id)
            job.failure_reason = failure_message
            raise ValueError(failure_message)


def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path,
                     downloader_job=job.id)
        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)
    except Exception:
        logger.exception("Exception caught while downloading batch.",
                         downloader_job=job.id)
        job.failure_reason = "Exception caught while downloading batch"
        raise
    finally:
        target_file.close()


def _extract_file(files: List[File], job: DownloaderJob) -> None:
    """Extract zip from temp directory and move to raw directory.

    Additionally this function sets the size_in_bytes field of each
    Batch in batches. To save database calls it does not save the
    batch itself since it will be saved soon when its status
    changes in utils.end_job.
    """
    # zip_path and local_dir should be common to all batches in the group
    job_dir = utils.JOB_DIR_PREFIX + str(job.id)
    zip_path = files[0].get_temp_download_path(job_dir)
    local_dir = files[0].get_temp_dir(job_dir)
    dirs_to_clean = set()

    logger.debug("Extracting %s", zip_path, downloader_job=job.id)

    try:
        zip_ref = zipfile.ZipFile(zip_path, "r")
        zip_ref.extractall(local_dir)

        for file in files:
            batch_directory = file.get_temp_dir(job_dir)
            raw_file_location = file.get_temp_pre_path(job_dir)

            # The platform is part of the batch's location so if the
            # batches in this job have different platforms then some
            # of them need to be moved to the directory corresponding
            # to thier platform.
            if local_dir != batch_directory:
                os.makedirs(batch_directory, exist_ok=True)
                dirs_to_clean.add(batch_directory)
                incorrect_location = os.path.join(local_dir, file.name)
                os.rename(incorrect_location, raw_file_location)

            file.size_in_bytes = os.path.getsize(raw_file_location)
            file.save()
            file.upload_raw_file(job_dir)
    except Exception:
        logger.exception("Exception caught while extracting %s",
                         zip_path,
                         downloader_job=job.id)
        job.failure_reason = "Exception caught while extracting " + zip_path
        raise
    finally:
        zip_ref.close()
        file.remove_temp_directory(job_dir)
        for directory in dirs_to_clean:
            shutil.rmtree(directory)


def download_array_express(job_id: int) -> None:
    """The main function for the Array Express Downloader.

    Downloads a single zip file containing the .PCL files representing
    samples relating to a single experiement stored in
    ArrayExpress. Each of these files is a separate Batch, so the file
    is unzipped and then each Batch's data is stored in Temporary
    Storage.
    """
    job = utils.start_job(job_id)
    batches = job.batches.all()
    success = True
    job_dir = utils.JOB_DIR_PREFIX + str(job_id)

    if batches.count() > 0:
        files = File.objects.filter(batch__in=batches)
        target_directory = files[0].get_temp_dir(job_dir)
        os.makedirs(target_directory, exist_ok=True)
        target_file_path = files[0].get_temp_download_path(job_dir)
        download_url = files[0].download_url
    else:
        logger.error("No batches found.",
                     downloader_job=job_id)
        success = False

    if success:
        try:
            _verify_batch_grouping(files, job)

            # The files for all of the batches in the grouping are
            # contained within the same zip file. Therefore only
            # download the one.
            _download_file(download_url, target_file_path, job)
            _extract_file(files, job)
        except Exception:
            # Exceptions are already logged and handled.
            # Just need to mark the job as failed.
            success = False

    if success:
        logger.debug("File %s downloaded and extracted successfully.",
                     download_url,
                     downloader_job=job_id)

    utils.end_job(job, batches, success)

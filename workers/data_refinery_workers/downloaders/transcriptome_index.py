from __future__ import absolute_import, unicode_literals

import os
import shutil
import urllib.request
from typing import List
from contextlib import closing

from data_refinery_common.models import (
    DownloaderJob,
    ProcessorJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJobOriginalFileAssociation,
)
from data_refinery_common.message_queue import send_job
from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.downloaders import utils

logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
CHUNK_SIZE = 1024 * 256 # chunk_size is in bytes


def _download_file(download_url: str, file_path: str, job: DownloaderJob) -> None:
    """Download the file via FTP.

    I spoke to Erin from Ensembl about ways to improve this. They're looking into it,
    but have decided against adding an Aspera endpoint.

    She suggested using `rsync`, we could try shelling out to that.

    """
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
    return True

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

    file_assocs = DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=job)
    files_to_process = []
    for assoc in file_assocs:
        original_file = assoc.original_file

        if original_file.is_archive:
            filename_species = ''.join(original_file.source_filename.split('.')[:-2])
        else:
            # Does this ever happen?
            filename_species = ''.join(original_file.source_filename.split('.')[:-1])

        os.makedirs(LOCAL_ROOT_DIR + '/' + filename_species, exist_ok=True)
        dl_file_path = LOCAL_ROOT_DIR + '/' + filename_species + '/' + original_file.source_filename
        success = _download_file(original_file.source_url, dl_file_path, job)

        if success:
            original_file.is_downloaded = True
            original_file.absolute_file_path = dl_file_path
            original_file.filename = original_file.source_filename
            original_file.is_archive = True
            original_file.has_raw = True
            original_file.calculate_size()
            original_file.calculate_sha1()
            original_file.save()
            files_to_process.append(original_file)
        else:
            logger.error("Problem during download",
                         url=original_file.source_url,
                         original_file_id=original_file.id,
                         downloader_job=job_id)

    if success:
        logger.debug("Files downloaded successfully.",
                     downloader_job=job_id)

    utils.end_downloader_job(job, success)
    create_long_and_short_processor_jobs(files_to_process)

def create_long_and_short_processor_jobs(files_to_process):
    """ Creates two processor jobs for the files needed for this transcriptome"""

    processor_job_long = ProcessorJob()
    processor_job_long.pipeline_applied = "TRANSCRIPTOME_INDEX_LONG"
    processor_job_long.ram_amount = 8192
    processor_job_long.save()

    for original_file in files_to_process:

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job_long
        assoc.save()

    send_job(ProcessorPipeline[processor_job_long.pipeline_applied], processor_job_long)

    processor_job_short = ProcessorJob()
    processor_job_short.pipeline_applied = "TRANSCRIPTOME_INDEX_SHORT"
    processor_job_long.ram_amount = 8192
    processor_job_short.save()

    for original_file in files_to_process:

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = original_file
        assoc.processor_job = processor_job_short
        assoc.save()

    send_job(ProcessorPipeline[processor_job_short.pipeline_applied], processor_job_short)


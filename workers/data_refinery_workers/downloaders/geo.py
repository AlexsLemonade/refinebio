import gzip
import os
import shutil
import subprocess
import tarfile
import time
import urllib.request
import re

from contextlib import closing
from typing import List, Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    Sample,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.downloaders import utils
from data_refinery_common.job_management import create_processor_jobs_for_original_files


logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# chunk_size is in bytes
CHUNK_SIZE = 1024 * 256


def _download_file(download_url: str, file_path: str, job: DownloaderJob, force_ftp=False) -> None:
    """ Download a file from GEO via Aspera unless `force_ftp` is True
    """

    # Ensure directory exists
    os.makedirs(file_path.rsplit('/', 1)[0], exist_ok=True)

    if not force_ftp:
        return _download_file_aspera(download_url=download_url, downloader_job=job, target_file_path=file_path)
    else:
        try:
            logger.debug("Downloading file from %s to %s.",
                         download_url,
                         file_path,
                         downloader_job=job.id)

            # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
            urllib.request.urlcleanup()

            target_file = open(file_path, "wb")
            with closing(urllib.request.urlopen(download_url)) as request:
                shutil.copyfileobj(request, target_file, CHUNK_SIZE)

            urllib.request.urlcleanup()
        except Exception:
            logger.exception("Exception caught while downloading file.",
                             downloader_job=job.id)
            job.failure_reason = "Exception caught while downloading file"
            raise
        finally:
            target_file.close()

        return True


def _download_file_aspera(download_url: str,
                          downloader_job: DownloaderJob,
                          target_file_path: str,
                          attempt=0) -> bool:
    """ Download a file to a location using Aspera by shelling out to the `ascp` client. """

    try:
        logger.debug("Downloading file from %s to %s via Aspera.",
                     download_url,
                     target_file_path,
                     downloader_job=downloader_job.id)

        ascp = ".aspera/cli/bin/ascp"
        key = ".aspera/cli/etc/asperaweb_id_dsa.openssh"
        url = download_url
        user = "anonftp"
        ftp = "ftp-trace.ncbi.nlm.nih.gov"
        if url.startswith("ftp://"):
            url = url.replace("ftp://", "")
        url = url.replace(ftp, "").replace('ftp.ncbi.nlm.nih.gov', '')

        # Resume level 1, use encryption, unlimited speed
        command_str = "{} -i {} -k1 -T {}@{}:{} {}".format(ascp, key, user, ftp, url, target_file_path)
        formatted_command = command_str.format(src=download_url,
                                               dest=target_file_path)
        completed_command = subprocess.run(formatted_command.split(),
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

        # Something went wrong! Else, just fall through to returning True.
        if completed_command.returncode != 0:

            stderr = completed_command.stderr.decode().strip()
            logger.debug("Shell call of `%s` to ascp failed with error message: %s",
                         formatted_command,
                         stderr,
                         downloader_job=downloader_job.id)

            # Sometimes, GEO fails mysteriously.
            # Wait a few minutes and try again.
            if attempt >= 5:
                downloader_job.failure_reason = stderr
                logger.error("All attempts to download accession via ascp failed: %s\nCommand was: %s",
                             stderr,
                             formatted_command,
                             downloader_job=downloader_job.id)
                return False
            else:
                time.sleep(30)
                return _download_file_aspera(download_url,
                                             downloader_job,
                                             target_file_path,
                                             attempt + 1
                                             )
    except Exception:
        logger.exception("Exception caught while downloading file from the URL via Aspera: %s",
                         download_url,
                         downloader_job=downloader_job.id)
        downloader_job.failure_reason = ("Exception caught while downloading "
                                         "file from the URL via Aspera: {}").format(download_url)
        return False

    # If Aspera has given a zero-byte file for some reason, let's back off and retry.
    if os.path.getsize(target_file_path) < 1:
        os.remove(target_file_path)
        if attempt > 5:
            downloader_job.failure_reason = "Got zero byte file from aspera after 5 attempts."
            return False

        logger.error("Got zero byte ascp download for target, retrying.",
                     target_url=download_url,
                     downloader_job=downloader_job.id)
        time.sleep(10)
        return _download_file_aspera(download_url,
                                     downloader_job,
                                     target_file_path,
                                     attempt + 1
                                     )
    return True

class FileExtractionError(Exception):
    def __init__(self, file_path, original_error):
        self.file_path = file_path
        self.original_error = original_error

class ArchivedFile:
    """ Contains utility functions to enumerate the files inside an archive that was downloaded from GEO. """
    def __init__(self, downloaded_file_path, parent_archive = None):
        self.file_path = downloaded_file_path
        self.parent_archive = parent_archive

        # thanks to https://stackoverflow.com/a/541394/763705
        self.filename = os.path.basename(self.file_path)
        self.extension = os.path.splitext(self.file_path)[1]

    def parent_is_archive(self):
        """ Returns true if this file came from an archive """
        return self.parent_archive and self.parent_archive.is_archive()

    def experiment_accession_code(self):
        """ Tries to get an experiment accession code from the file name.
        GEO experiment accession codes start with GSE and have between 5 and 9 characters """
        match = re.match(r'(GSE\d{2,6})', self.filename)
        if match: return match.group(0)
        return None

    def sample_accession_code(self):
        """ Tries to get a sample accession code from the file name.
        GEO sample accession codes start with GSM and have between 6 and 10 characters """
        match = re.match(r'(GSM\d{4,7})', self.filename)
        if match: return match.group(0)
        return None

    def get_sample(self):
        """ Tries to find the sample associated with this file, and returns None if unable. """
        return Sample.objects.filter(accession_code=self.sample_accession_code()).first()

    def is_processable(self):
        """ There're some known file patterns that are found in GEO that we know we can ignore. """
        if re.match(r'(GSE\d{2,6})_family.xml', self.filename):
            return False

        # ignore platform files
        if re.match(r'(GPL\d{3,5})', self.filename):
            return False

        return True

    def is_archive(self):
        return self.extension.lower() in ['.tar', '.tgz', '.gz']

    def get_files(self):
        if not self.is_archive():
            yield self
        else:
            # for archives extract them and enumerate all the files inside
            for path in self._extract_files():
                archived_file = ArchivedFile(path, self)
                for file in archived_file.get_files():
                    yield file

    def _extract_files(self):
        logger.debug("Extracting %s!", self.file_path, file_path=self.file_path)

        try:
            if '.tar' == self.extension:
                return self._extract_tar()
            elif '.tgz' == self.extension:
                return self._extract_tgz()
            elif '.gz' == self.extension:
                return self._extract_gz()
        except Exception as e:
            logger.exception("While extracting %s caught exception %s", self.file_path, str(e), file_path=self.file_path)
            raise FileExtractionError(self.file_path, e)
        
        raise FileExtractionError(self.file_path, 'Unknown archive file format.')

    def _get_absolute_path(self):
        """ This returns the path where this file would be extracted if it's an archive """
        if self.experiment_accession_code():
            return LOCAL_ROOT_DIR + '/' + self.experiment_accession_code() + '/raw/'
        elif self.sample_accession_code():
            return LOCAL_ROOT_DIR + '/' + self.sample_accession_code() + '/raw/'
        else:
            return LOCAL_ROOT_DIR + '/' + self.filename + '/raw/'

    def _extract_tar(self) -> List[str]:
        """ Extract tar and return a list of the raw files. """
        # This is technically an unsafe operation.
        # However, we're trusting GEO as a data source.
        abs_with_code_raw = self._get_absolute_path()

        with tarfile.TarFile(self.file_path, "r") as zip_ref:
            zip_ref.extractall(abs_with_code_raw)
            extracted_files = zip_ref.getnames() # https://docs.python.org/3/library/tarfile.html#tarfile.TarFile.getnames

        return [(abs_with_code_raw + f) for f in extracted_files]

    def _extract_tgz(self) -> List[str]:
        """Extract tgz and return a list of the raw files."""
        abs_with_code_raw = self._get_absolute_path()
        
        extracted_filepath = self.file_path.replace('.tgz', '.tar')

        with gzip.open(self.file_path, 'rb') as f_in:
            with open(extracted_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        with tarfile.TarFile(extracted_filepath, "r") as zip_ref:
            zip_ref.extractall(abs_with_code_raw)
            extracted_files = zip_ref.getnames() # https://docs.python.org/3/library/tarfile.html#tarfile.TarFile.getnames

        return [(abs_with_code_raw + f) for f in extracted_files]

    def _extract_gz(self) -> List[str]:
        """Extract gz and return a list of the raw files."""
        extracted_filepath = self.file_path.replace('.gz', '')
        with gzip.open(self.file_path, 'rb') as f_in:
            with open(extracted_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return [extracted_filepath]



def download_geo(job_id: int) -> None:
    """The main function for the GEO Downloader.

    Downloads a single tar file containing the files representing
    samples relating to a single experiment stored in
    GEO.
    """
    job = utils.start_job(job_id)
    accession_code = job.accession_code
    original_file = job.original_files.first()

    if not original_file:
        job.failure_reason = "No files associated with the job."
        logger.error("No files associated with the job.", downloader_job=job_id)
        utils.end_downloader_job(job, success=False)
        return

    url = original_file.source_url
    related_samples = original_file.samples.exclude(technology='RNA-SEQ')

    # First, download the sample archive URL.
    # Then, unpack all the ones downloaded.
    # Then create processor jobs!

    # The files for all of the samples are contained within the same zip file. Therefore only
    # download the one.
    os.makedirs(LOCAL_ROOT_DIR + '/' + accession_code, exist_ok=True)
    dl_file_path = LOCAL_ROOT_DIR + '/' + accession_code + '/' + url.split('/')[-1]

    logger.debug("Starting to download: " + url, job_id=job_id, accession_code=accession_code)
    _download_file(url, dl_file_path, job)
    original_file.absolute_file_path = dl_file_path
    original_file.is_downloaded = True
    original_file.save()

    has_raw = True
    unpacked_sample_files = []

    try:
        # enumerate all files inside the archive
        archived_files = list(ArchivedFile(dl_file_path).get_files())
    except FileExtractionError as e:
        job.failure_reason = e
        logger.exception("Error occurred while extracting file.", path=dl_file_path, exception=str(e))
        utils.end_downloader_job(job, success=False)
        return

    for og_file in archived_files:
        sample = og_file.get_sample()

        # We don't want RNA-Seq data from GEO:
        # https://github.com/AlexsLemonade/refinebio/issues/966
        if sample and sample.technology == 'RNA-SEQ':
            logger.warn("RNA-Seq sample found in GEO downloader job.", sample=sample)
            continue

        if not sample and (not og_file.is_processable() or og_file.experiment_accession_code() != accession_code):
            # skip the files that we know are not processable and can't be associated with a sample
            # also skip the files were we couldn't find a sample and they don't mention the current experiment            
            continue

        potential_existing_file = OriginalFile.objects.filter(
            source_filename=original_file.source_filename,
            filename=og_file.filename,
            is_archive=False
        ).first()
        if potential_existing_file:
            # We've already created this record, let's see if we actually
            # needed to download it or if we just got it because we needed
            # a file in the same archive.
            if potential_existing_file.needs_processing():
                if not potential_existing_file.is_downloaded:
                    potential_existing_file.is_downloaded = True
                    potential_existing_file.save()

                unpacked_sample_files.append(potential_existing_file)
            continue

        # Then this is a new file and we should create an original file for it
        actual_file = OriginalFile()
        actual_file.is_downloaded = True
        actual_file.is_archive = False
        actual_file.absolute_file_path = og_file.file_path
        actual_file.filename = og_file.filename
        actual_file.calculate_size()
        actual_file.calculate_sha1()
        actual_file.has_raw = True
        actual_file.source_url = original_file.source_url
        actual_file.source_filename = original_file.source_filename
        actual_file.save()

        # try to see if the file should be associated with a sample
        if sample:
            original_file_sample_association = OriginalFileSampleAssociation()
            original_file_sample_association.sample = sample
            original_file_sample_association.original_file = actual_file
            original_file_sample_association.save()
        else:
            # if not, we can associate this file with all samples in the experiment
            for sample in related_samples:
                original_file_sample_association = OriginalFileSampleAssociation()
                original_file_sample_association.sample = sample
                original_file_sample_association.original_file = actual_file
                original_file_sample_association.save()

        unpacked_sample_files.append(actual_file)

    if len(unpacked_sample_files) > 0:
        success = True
        logger.debug("File downloaded and extracted successfully.",
                     url=url,
                     dl_file_path=dl_file_path,
                     downloader_job=job_id)
    else:
        success = False
        logger.info("Unable to extract any files.",
                    url=url,
                    dl_file_path=dl_file_path,
                    downloader_job=job_id)
        job.failure_reason = "Failed to extract any downloaded files."

    if success:
        create_processor_jobs_for_original_files(unpacked_sample_files, job)

    if original_file.is_archive:
        original_file.delete_local_file()

    utils.end_downloader_job(job, success)

    return success

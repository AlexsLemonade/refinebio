from typing import Dict, Set, List
from urllib.parse import urlparse
import csv
import nomad
import io
import os
import re
import requests

from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from retrying import retry

import hashlib
from functools import partial

# Found: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
METADATA_URL = "http://169.254.169.254/latest/meta-data"
INSTANCE_ID = None
SUPPORTED_MICROARRAY_PLATFORMS = None
SUPPORTED_RNASEQ_PLATFORMS = None
READABLE_PLATFORM_NAMES = None

def get_env_variable(var_name: str, default: str=None) -> str:
    """ Get an environment variable or return a default value """
    try:
        return os.environ[var_name]
    except KeyError:
        if default:
            return default
        error_msg = "Set the %s environment variable" % var_name
        raise ImproperlyConfigured(error_msg)


def get_env_variable_gracefully(var_name: str, default: str=None) -> str:
    """
    Get an environment variable, or return a default value, but always fail gracefully and return
    something rather than raising an ImproperlyConfigured error.
    """
    try:
        return os.environ[var_name]
    except KeyError:
        return default


def get_instance_id() -> str:
    """Returns the AWS instance id where this is running or "local"."""
    global INSTANCE_ID
    if INSTANCE_ID is None:
        if settings.RUNNING_IN_CLOUD:
            @retry(stop_max_attempt_number=3)
            def retrieve_instance_id():
                return requests.get(os.path.join(METADATA_URL, "instance-id")).text

            INSTANCE_ID = retrieve_instance_id()
        else:
            INSTANCE_ID = "local"

    return INSTANCE_ID


def get_worker_id() -> str:
    """Returns <instance_id>/<thread_id>."""
    return get_instance_id() + "/" + current_process().name


def get_volume_index(path='/home/user/data_store/VOLUME_INDEX') -> str:
    """ Reads the contents of the VOLUME_INDEX file, else returns default """

    if settings.RUNNING_IN_CLOUD:
        default = "-1"
    else:
        default = "0"

    try:
        with open(path, 'r') as f:
            v_id = f.read().strip()
            return v_id
    except Exception as e:
        # Our configured logger needs util, so we use the standard logging library for just this.
        import logging
        logger = logging.getLogger(__name__)
        logger.info(str(e))
        logger.info("Could not read volume index file, using default: " + str(default))

    return default

def get_nomad_jobs() -> list:
    """Calls nomad service and return all jobs"""
    try:
        nomad_host = get_env_variable("NOMAD_HOST")
        nomad_port = get_env_variable("NOMAD_PORT", "4646")
        nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)
        return nomad_client.jobs.get_jobs()
    except nomad.api.exceptions.BaseNomadException:
        # Nomad is not available right now
        return []

def get_active_volumes() -> Set[str]:
    """Returns a Set of indices for volumes that are currently mounted.

    These can be used to determine which jobs would actually be able
    to be placed if they were queued up.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

    volumes = set()
    try:
        for node in nomad_client.nodes.get_nodes():
            node_detail = nomad_client.node.get_node(node["ID"])
            if 'Status' in node_detail and node_detail['Status'] == 'ready' \
               and 'Meta' in node_detail and 'volume_index' in node_detail['Meta']:
                volumes.add(node_detail['Meta']['volume_index'])
    except nomad.api.exceptions.BaseNomadException:
        # Nomad is down, return the empty set.
        pass

    return volumes


def get_supported_microarray_platforms(platforms_csv: str="config/supported_microarray_platforms.csv"
                                       ) -> list:
    """
    Loads our supported microarray platforms file and returns a list of dictionaries
    containing the internal accession, the external accession, and a boolean indicating
    whether or not the platform supports brainarray.
    CSV must be in the format:
    Internal Accession | External Accession | Supports Brainarray
    """
    global SUPPORTED_MICROARRAY_PLATFORMS
    if SUPPORTED_MICROARRAY_PLATFORMS is not None:
        return SUPPORTED_MICROARRAY_PLATFORMS

    SUPPORTED_MICROARRAY_PLATFORMS = []
    with open(platforms_csv) as platforms_file:
        reader = csv.reader(platforms_file)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            external_accession = line[1]
            is_brainarray = True if line[2] == 'y' else False
            SUPPORTED_MICROARRAY_PLATFORMS.append({"platform_accession": line[0],
                                                   "external_accession": external_accession,
                                                   "is_brainarray": is_brainarray})

            # A-GEOD-13158 is the same platform as GPL13158 and this
            # pattern is generalizable. Since we don't want to have to
            # list a lot of platforms twice just with different prefixes,
            # we just convert them and add them to the list.
            if external_accession[:6] == "A-GEOD":
                converted_accession = external_accession.replace("A-GEOD-", "GPL")
                SUPPORTED_MICROARRAY_PLATFORMS.append({"platform_accession": line[0],
                                                       "external_accession": converted_accession,
                                                       "is_brainarray": is_brainarray})

            # Our list of supported platforms contains both A-GEOD-*
            # and GPL*, so convert both ways.
            if external_accession[:3] == "GPL":
                converted_accession = external_accession.replace("GPL", "A-GEOD-")
                SUPPORTED_MICROARRAY_PLATFORMS.append({"platform_accession": line[0],
                                                       "external_accession": converted_accession,
                                                       "is_brainarray": is_brainarray})

    return SUPPORTED_MICROARRAY_PLATFORMS


def get_supported_rnaseq_platforms(platforms_list: str="config/supported_rnaseq_platforms.txt"
                                   ) -> list:
    """
    Returns a list of RNASeq platforms which are currently supported.
    """
    global SUPPORTED_RNASEQ_PLATFORMS
    if SUPPORTED_RNASEQ_PLATFORMS is not None:
        return SUPPORTED_RNASEQ_PLATFORMS

    SUPPORTED_RNASEQ_PLATFORMS = []
    with open(platforms_list) as platforms_file:
        for line in platforms_file:
            SUPPORTED_RNASEQ_PLATFORMS.append(line.strip())

    return SUPPORTED_RNASEQ_PLATFORMS


def get_readable_affymetrix_names(mapping_csv: str="config/readable_affymetrix_names.csv") -> Dict:
    """
    Loads the mapping from human readble names to internal accessions for Affymetrix platforms.
    CSV must be in the format:
    Readable Name | Internal Accession
    Returns a dictionary mapping from internal accessions to human readable names.
    """
    global READABLE_PLATFORM_NAMES

    if READABLE_PLATFORM_NAMES is not None:
        return READABLE_PLATFORM_NAMES

    READABLE_PLATFORM_NAMES = {}
    with open(mapping_csv, encoding='utf-8') as mapping_file:
        reader = csv.reader(mapping_file, )
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            READABLE_PLATFORM_NAMES[line[1]] = line[0]

    return READABLE_PLATFORM_NAMES


def get_internal_microarray_accession(accession_code):
    platforms = get_supported_microarray_platforms()

    for platform in platforms:
        if platform['external_accession'] == accession_code:
            return platform['platform_accession']
        elif platform['platform_accession'] == accession_code:
            return platform['platform_accession']

    return None

def get_normalized_platform(external_accession):
    """
    Handles a weirdo cases, where external_accessions in the format
        hugene10stv1 -> hugene10st
    """

    matches = re.findall(r"stv\d$", external_accession)
    for match in matches:
        external_accession = external_accession.replace(match, 'st')

    return external_accession

def parse_s3_url(url):
    """
    Parses S3 URL.
    Returns bucket (domain) and file (full path).
    """
    bucket = ''
    path = ''
    if url:
        result = urlparse(url)
        bucket = result.netloc
        path = result.path.strip('/')
    return bucket, path

def get_s3_url(s3_bucket: str, s3_key: str) -> str:
    """
    Calculates the s3 URL for a file from the bucket name and the file key.
    """
    return "%s.s3.amazonaws.com/%s" % (s3_bucket, s3_key)

def calculate_file_size(absolute_file_path):
    return os.path.getsize(absolute_file_path)

def calculate_sha1(absolute_file_path):
    hash_object = hashlib.sha1()
    with open(absolute_file_path, mode='rb') as open_file:
        for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b''):
            hash_object.update(buf)

    return hash_object.hexdigest()

def get_sra_download_url(run_accession, protocol="fasp"):
    """Try getting the sra-download URL from CGI endpoint"""
    #Ex: curl --data "acc=SRR6718414&accept-proto=fasp&version=2.0" https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi
    cgi_url = "https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi"
    data = "acc=" + run_accession + "&accept-proto=" + protocol + "&version=2.0"
    try:
        resp = requests.post(cgi_url, data=data)
    except Exception as e:
        # Our configured logger needs util, so we use the standard logging library for just this.
        import logging
        logger = logging.getLogger(__name__)
        logger.exception("Bad CGI request!: " + str(cgi_url) + ", " + str(data))
        return None

    if resp.status_code != 200:
        # This isn't on the new servers
        return None
    else:
        try:
            # From: '#2.0\nsrapub|DRR002116|2324796808|2013-07-03T05:51:55Z|50964cfc69091cdbf92ea58aaaf0ac1c||fasp://dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/dra0/DRR/000002/DRR002116|200|ok\n'
            # To:  'dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/dra0/DRR/000002/DRR002116'

            # Sometimes, the responses from names.cgi makes no sense at all on a per-accession-code basis. This helps us handle that.
            # $ curl --data "acc=SRR5818019&accept-proto=fasp&version=2.0" https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi
            # 2.0\nremote|SRR5818019|434259775|2017-07-11T21:32:08Z|a4bfc16dbab1d4f729c4552e3c9519d1|||400|Only 'https' protocol is allowed for this object
            protocol_header = protocol + '://'
            sra_url = resp.text.split('\n')[1].split('|')[6]
            return sra_url
        except Exception as e:
            # Our configured logger needs util, so we use the standard logging library for just this.
            import logging
            logger = logging.getLogger(__name__)
            logger.exception("Error parsing CGI response: " + str(cgi_url) + " " + str(data) + " " + str(resp.text))
            return None


def get_fasp_sra_download(run_accession: str):
    """Get an URL for SRA using the FASP protocol.

    These URLs should not actually include the protcol."""
    full_url = get_sra_download_url(run_accession, 'fasp')
    if full_url:
        sra_url = full_url.split('fasp://')[1]
        return sra_url
    else:
        return None

def get_https_sra_download(run_accession: str):
    """Get an HTTPS URL for SRA."""
    return get_sra_download_url(run_accession, 'https')

def load_blacklist(blacklist_csv: str="config/RNASeqRunBlackList.csv"):
    """ Loads the SRA run blacklist """

    blacklisted_samples = []
    with open(blacklist_csv, encoding='utf-8') as blacklist_file:
        reader = csv.reader(blacklist_file, )
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            blacklisted_samples.append(line[0].strip())

    return blacklisted_samples


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

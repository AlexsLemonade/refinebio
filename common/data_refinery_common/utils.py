import csv
import os
import re
import requests
from urllib.parse import urlparse
from typing import Dict

from billiard import current_process
from django.core.exceptions import ImproperlyConfigured
from retrying import retry

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
        if get_env_variable("RUNNING_IN_CLOUD") == "True":
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


def get_volume_index(default="0", path='/home/user/data_store/VOLUME_INDEX') -> str:
    """ Reads the contents of the VOLUME_INDEX file, else returns default """

    try:
        with open(path, 'r') as f:
            v_id = f.read().strip()
            return v_id
    except Exception as e:
        # Logger needs util, so we do this at runtime
        from data_refinery_common.logging import get_and_configure_logger
        logger = get_and_configure_logger(__name__)
        logger.info("Could not read volume index file, using default")

    return default


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

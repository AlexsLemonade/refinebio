from typing import Dict, Set
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

def get_fasp_sra_download(run_accession: str):
    """Try getting the sra-download URL from CGI endpoint"""
    #Ex: curl --data "acc=SRR6718414&accept-proto=fasp&version=2.0" https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi
    cgi_url = "https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi"
    data = "acc=" + run_accession + "&accept-proto=fasp&version=2.0"
    try:
        resp = requests.post(cgi_url, data=data)
    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.exception("Bad FASP CGI request", data=data)
        return None

    if resp.status_code != 200:
        # This isn't on the new FASP servers
        return None
    else:
        try:
            # From: '#2.0\nsrapub|DRR002116|2324796808|2013-07-03T05:51:55Z|50964cfc69091cdbf92ea58aaaf0ac1c||fasp://dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/dra0/DRR/000002/DRR002116|200|ok\n'
            # To:  'dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/dra0/DRR/000002/DRR002116'
            sra_url = resp.text.split('\n')[1].split('|')[6].split('fasp://')[1]
            return sra_url
        except Exception as e:
            # Our configured logger needs util, so we use the standard logging library for just this.
            import logging
            logger = logging.getLogger(__name__)
            logger.exception("Bad FASP CGI response", data=data, text=resp.text)
            return None

def has_original_file_been_processed(original_file) -> bool:
    """Returns True if original_file is from SRA and has been processed

    Returns False otherwise.

    original_file must have source_database == SRA, otherwise an
    exception will be raised. This is because SRA does not have more
    than one sample per original file, which is a precondition for
    this function being correct. The reason for this is that otherwise
    an archive file could have many samples which are unprocessed, but
    if there's one then we would return True here which would be
    misleading.
    """
    sample = original_file.samples().first()

    if not sample:
        return False

    if sample.source_database != "SRA":
        raise Exception(("has_original_file_been_processed called on an OriginalFile that"
                         " was not from SRA! This is unsupported behavior!"))

    for computed_file in sample.computed_files.all():
            if compted_file.s3_bucket and compted_file.s3_key:
                return True

    return False

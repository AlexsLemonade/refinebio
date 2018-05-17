import csv
import os
import requests
from urllib.parse import urlparse
from typing import Dict

from billiard import current_process
from django.core.exceptions import ImproperlyConfigured
from retrying import retry

# Found: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
METADATA_URL = "http://169.254.169.254/latest/meta-data"
INSTANCE_ID = None


def get_env_variable(var_name: str, default:str=None) -> str:
    """ Get an environment variable or return a default value """
    try:
        return os.environ[var_name]
    except KeyError:
        if default:
            return default
        error_msg = "Set the %s environment variable" % var_name
        raise ImproperlyConfigured(error_msg)


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

def get_supported_microarray_platforms(platforms_csv:str="supported_microarray_platforms.csv") -> list:
    """
    Loads our supported microarray platforms file and returns a list of dictionaries
    containing the internal accession, the external accession, and a boolean indicating
    whether or not the platform supports brainarray.
    CSV must be in the format:
    Internal Accession | External Accession | Supports Brainarray
    """
    supported_platforms = []
    with open(platforms_csv) as platforms_file:
        reader = csv.reader(platforms_file)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            supported_platforms.append({"platform_accession": line[0],
                                        "external_accession": line[1],
                                        "is_brainarray": True if line[2] == 'y' else False})

    return supported_platforms

def get_supported_rnaseq_platforms(platforms_list:str="supported_rnaseq_platforms.txt") -> list:
    """
    Returns a list of RNASeq platforms which are currently supported.
    """
    supported_platforms = []
    with open(platforms_list) as platforms_file:
        for line in platforms_file:
            supported_platforms.append(line.strip())

    return supported_platforms

def get_readable_platform_names(mapping_csv:str="readable_platform_names.csv") -> Dict:
    """
    Loads the mapping from human readble names to internal accessions for Microarray platforms.
    CSV must be in the format:
    Internal Accession | External Accession | Supports Brainarray
    Returns a dictionary mapping from internal accessions to human readable names.
    """
    names_mapping = {}
    with open(mapping_csv) as mapping_file:
        reader = csv.reader(mapping_file)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            names_mapping[line[1]] = line[0]

    return names_mapping

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


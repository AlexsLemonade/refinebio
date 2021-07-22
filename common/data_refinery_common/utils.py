import csv
import hashlib
import io
import os
import random
import re
import time
from functools import partial
from multiprocessing import current_process
from typing import Dict
from urllib.parse import urlparse

from django.conf import settings
from django.core.exceptions import ImproperlyConfigured

import requests
from retrying import retry

from data_refinery_common.performant_pagination.pagination import PerformantPaginator

# Found: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
METADATA_URL = "http://169.254.169.254/latest/meta-data"
INSTANCE_ID = None
SUPPORTED_MICROARRAY_PLATFORMS = None
SUPPORTED_RNASEQ_PLATFORMS = None
READABLE_PLATFORM_NAMES = None
ANNOTATION_PACKAGE_OVERRIDES = None


def get_env_variable(var_name: str, default: str = None) -> str:
    """ Get an environment variable or return a default value """
    try:
        return os.environ[var_name]
    except KeyError:
        if default:
            return default
        error_msg = "Set the %s environment variable" % var_name
        raise ImproperlyConfigured(error_msg)


def get_env_variable_gracefully(var_name: str, default: str = None) -> str:
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


def get_supported_microarray_platforms(
    platforms_csv: str = "config/supported_microarray_platforms.csv",
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
            if reader.line_num == 1:
                continue

            external_accession = line[1]
            is_brainarray = True if line[2] == "y" else False
            SUPPORTED_MICROARRAY_PLATFORMS.append(
                {
                    "platform_accession": line[0],
                    "external_accession": external_accession,
                    "is_brainarray": is_brainarray,
                }
            )

            # A-GEOD-13158 is the same platform as GPL13158 and this
            # pattern is generalizable. Since we don't want to have to
            # list a lot of platforms twice just with different prefixes,
            # we just convert them and add them to the list.
            if external_accession[:6] == "A-GEOD":
                converted_accession = external_accession.replace("A-GEOD-", "GPL")
                SUPPORTED_MICROARRAY_PLATFORMS.append(
                    {
                        "platform_accession": line[0],
                        "external_accession": converted_accession,
                        "is_brainarray": is_brainarray,
                    }
                )

            # Our list of supported platforms contains both A-GEOD-*
            # and GPL*, so convert both ways.
            if external_accession[:3] == "GPL":
                converted_accession = external_accession.replace("GPL", "A-GEOD-")
                SUPPORTED_MICROARRAY_PLATFORMS.append(
                    {
                        "platform_accession": line[0],
                        "external_accession": converted_accession,
                        "is_brainarray": is_brainarray,
                    }
                )

    return SUPPORTED_MICROARRAY_PLATFORMS


def get_supported_rnaseq_platforms(
    platforms_list: str = "config/supported_rnaseq_platforms.txt",
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


def get_readable_affymetrix_names(
    mapping_csv: str = "config/readable_affymetrix_names.csv",
) -> Dict:
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
    with open(mapping_csv, encoding="utf-8") as mapping_file:
        reader = csv.reader(mapping_file,)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num == 1:
                continue

            READABLE_PLATFORM_NAMES[line[1]] = line[0]

    return READABLE_PLATFORM_NAMES


def get_affymetrix_annotation_package_name_overrides(
    overrides_csv: str = "config/affymetrix_annotation_package_name_overrides.csv",
) -> Dict:
    """
    Loads the mapping from annotation package name to internal accession for Affymetrix platforms.
    CSV must be in the format:
    Annotation Package Name | Inernal Accession
    """
    global ANNOTATION_PACKAGE_OVERRIDES

    if ANNOTATION_PACKAGE_OVERRIDES is not None:
        return ANNOTATION_PACKAGE_OVERRIDES

    ANNOTATION_PACKAGE_OVERRIDES = {}
    with open(overrides_csv, encoding="utf-8") as overrides_file:
        reader = csv.reader(overrides_file,)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num == 1:
                continue

            ANNOTATION_PACKAGE_OVERRIDES[line[1]] = line[0]

    return ANNOTATION_PACKAGE_OVERRIDES


def get_internal_microarray_accession(accession_code):
    platforms = get_supported_microarray_platforms()

    for platform in platforms:
        if platform["external_accession"] == accession_code:
            return platform["platform_accession"]
        elif platform["platform_accession"] == accession_code:
            return platform["platform_accession"]

    return None


def get_normalized_platform(external_accession):
    """
    Handles a weirdo cases, where external_accessions in the format
        hugene10stv1 -> hugene10st
    """

    matches = re.findall(r"stv\d$", external_accession)
    for match in matches:
        external_accession = external_accession.replace(match, "st")

    return external_accession


def parse_s3_url(url):
    """
    Parses S3 URL.
    Returns bucket (domain) and file (full path).
    """
    bucket = ""
    path = ""
    if url:
        result = urlparse(url)
        bucket = result.netloc
        path = result.path.strip("/")
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
    with open(absolute_file_path, mode="rb") as open_file:
        for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b""):
            hash_object.update(buf)


def calculate_md5(absolute_file_path):
    hash_object = hashlib.md5()
    with open(absolute_file_path, mode="rb") as open_file:
        for buf in iter(partial(open_file.read, io.DEFAULT_BUFFER_SIZE), b""):
            hash_object.update(buf)

    return hash_object.hexdigest()


def get_sra_download_url(run_accession, protocol="fasp"):
    """Try getting the sra-download URL from CGI endpoint"""
    # Ex: curl --data "acc=SRR6718414&accept-proto=fasp&version=2.0" \
    #   https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi
    cgi_url = "https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi"
    data = "acc=" + run_accession + "&accept-proto=" + protocol + "&version=2.0"
    try:
        resp = requests.post(cgi_url, data=data)
    except Exception:
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
            sra_url = resp.text.split("\n")[1].split("|")[6]
            return sra_url
        except Exception:
            # Our configured logger needs util
            # so we use the standard logging library for just this.
            import logging

            logger = logging.getLogger(__name__)
            exception_template = "Error parsing CGI response: {0} {1} {2}"
            logger.exception(exception_template.format(str(cgi_url), str(data), str(resp.text)))
            return None


def get_fasp_sra_download(run_accession: str):
    """Get an URL for SRA using the FASP protocol.

    These URLs should not actually include the protcol."""
    full_url = get_sra_download_url(run_accession, "fasp")
    if full_url:
        sra_url = full_url.split("fasp://")[1]
        return sra_url
    else:
        return None


def get_https_sra_download(run_accession: str):
    """Get an HTTPS URL for SRA."""
    return get_sra_download_url(run_accession, "https")


def load_blacklist(blacklist_csv: str = "config/RNASeqRunBlackList.csv"):
    """ Loads the SRA run blacklist """

    blacklisted_samples = []
    with open(blacklist_csv, encoding="utf-8") as blacklist_file:
        reader = csv.reader(blacklist_file,)
        for line in reader:
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num == 1:
                continue

            blacklisted_samples.append(line[0].strip())

    return blacklisted_samples


def queryset_page_iterator(queryset, page_size=2000):
    """ use the performant paginator to iterate over each page in a queryset """
    paginator = PerformantPaginator(queryset, page_size)
    page = paginator.page()
    while True:
        yield page.object_list

        if not page.has_next():
            break
        else:
            page = paginator.page(page.next_page_number())


def queryset_iterator(queryset, page_size=2000):
    """ use the performant paginator to iterate over a queryset """
    for page in queryset_page_iterator(queryset, page_size):
        for item in page:
            yield item


def download_file(
    download_url,
    target_file_path,
    retry=1,
    *,
    backoff_factor=8,
    max_retries=10,
    max_sleep_timeout=120  # 2 mins
):
    """
    Downloads the given url into `target_file_path`.

    The download will be retried `max_retries` times if it fails for any reason.
    We use the exponential backoff https://en.wikipedia.org/wiki/Exponential_backoff
    algorithm to add some randomness between retries.
    """
    try:
        # thanks to https://stackoverflow.com/a/39217788/763705
        with requests.get(download_url, stream=True) as r:
            with open(target_file_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=(4 * 1024 * 1024)):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)
    except Exception as e:
        if retry < max_retries:
            # After the retry-th failed attempt, retry downloading after
            # k*backoff_factor, where k is a random integer between 0 and 2^retry − 1.
            k = random.randint(0, 2 ** retry - 1)
            sleep_timeout = min(k * backoff_factor, max_sleep_timeout)
            time.sleep(sleep_timeout)

            # and retry downloading again
            download_file(download_url, target_file_path, retry + 1)
        else:
            raise e


class FileUtils:
    @staticmethod
    def is_archive(file_path):
        extension = FileUtils.get_extension(file_path)
        if not extension:
            return False

        return extension in [".tar", ".tgz", ".gz", ".zip"]

    @staticmethod
    def get_filename(file_path):
        return os.path.basename(file_path)

    @staticmethod
    def get_extension(file_path):
        if not file_path:
            return None

        return os.path.splitext(file_path)[1].lower()

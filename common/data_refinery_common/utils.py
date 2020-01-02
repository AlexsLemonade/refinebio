import csv
import hashlib
import io
import os
import re
from functools import partial
from itertools import groupby
from typing import Dict, Set
from urllib.parse import urlparse
import nomad
import requests
from django.conf import settings
from django.core.exceptions import ImproperlyConfigured
from retrying import retry

from data_refinery_common.performant_pagination.pagination import PerformantPaginator

# Found: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
METADATA_URL = "http://169.254.169.254/latest/meta-data"
INSTANCE_ID = None
SUPPORTED_MICROARRAY_PLATFORMS = None
SUPPORTED_RNASEQ_PLATFORMS = None
READABLE_PLATFORM_NAMES = None


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


def get_volume_index(path="/home/user/data_store/VOLUME_INDEX") -> str:
    """ Reads the contents of the VOLUME_INDEX file, else returns default """

    if settings.RUNNING_IN_CLOUD:
        default = "-1"
    else:
        default = "0"

    try:
        with open(path, "r") as f:
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
            if (
                "Status" in node_detail
                and node_detail["Status"] == "ready"
                and "Meta" in node_detail
                and "volume_index" in node_detail["Meta"]
            ):
                volumes.add(node_detail["Meta"]["volume_index"])
    except nomad.api.exceptions.BaseNomadException:
        # Nomad is down, return the empty set.
        pass

    return volumes


def get_active_volumes_detailed() -> Dict:
    """Returns the instance type and number of allocations (jobs) for each active volume

    These can be used to determine which jobs would actually be able
    to be placed if they were queued up.
    """
    nomad_host = get_env_variable("NOMAD_HOST")
    nomad_port = get_env_variable("NOMAD_PORT", "4646")
    nomad_client = nomad.Nomad(nomad_host, port=int(nomad_port), timeout=30)

    volumes = dict()
    try:
        for node in nomad_client.nodes.get_nodes():
            node_detail = nomad_client.node.get_node(node["ID"])
            allocations = len(nomad_client.node.get_allocations(node["ID"]))
            if (
                "Status" in node_detail
                and node_detail["Status"] == "ready"
                and "Meta" in node_detail
                and "volume_index" in node_detail["Meta"]
            ):
                volume_info = dict()
                volume_info["type"] = node_detail["Attributes"].get(
                    "platform.aws.instance-type", None
                )
                volume_info["allocations"] = allocations

                volumes[node_detail["Meta"]["volume_index"]] = volume_info
    except nomad.api.exceptions.BaseNomadException:
        # Nomad is down, return the empty dict.
        pass

    return volumes


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
            if reader.line_num is 1:
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
            if reader.line_num is 1:
                continue

            READABLE_PLATFORM_NAMES[line[1]] = line[0]

    return READABLE_PLATFORM_NAMES


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

    return hash_object.hexdigest()


def get_sra_download_url(run_accession, protocol="fasp"):
    """Try getting the sra-download URL from CGI endpoint"""
    # Ex: curl --data "acc=SRR6718414&accept-proto=fasp&version=2.0" https://www.ncbi.nlm.nih.gov/Traces/names/names.cgi
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
            protocol_header = protocol + "://"
            sra_url = resp.text.split("\n")[1].split("|")[6]
            return sra_url
        except Exception as e:
            # Our configured logger needs util, so we use the standard logging library for just this.
            import logging

            logger = logging.getLogger(__name__)
            logger.exception(
                "Error parsing CGI response: "
                + str(cgi_url)
                + " "
                + str(data)
                + " "
                + str(resp.text)
            )
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
            if reader.line_num is 1:
                continue

            blacklisted_samples.append(line[0].strip())

    return blacklisted_samples


def get_nomad_jobs_breakdown():
    jobs = get_nomad_jobs()
    parameterized_jobs = [job for job in jobs if job["ParameterizedJob"]]

    def get_job_type(job):
        return get_job_details(job)[0]

    def get_job_volume(job):
        return get_job_details(job)[1]

    # groupby must be executed on a sorted iterable https://docs.python.org/2/library/itertools.html#itertools.groupby
    sorted_jobs_by_type = sorted(filter(get_job_type, parameterized_jobs), key=get_job_type)
    aggregated_jobs_by_type = groupby(sorted_jobs_by_type, get_job_type)
    nomad_pending_jobs_by_type, nomad_running_jobs_by_type = _aggregate_nomad_jobs(
        aggregated_jobs_by_type
    )

    # To get the total jobs for running and pending, the easiest
    # AND the most efficient way is to sum up the stats we've
    # already partially summed up.
    nomad_running_jobs = sum(num_jobs for job_type, num_jobs in nomad_running_jobs_by_type.items())
    nomad_pending_jobs = sum(num_jobs for job_type, num_jobs in nomad_pending_jobs_by_type.items())

    sorted_jobs_by_volume = sorted(filter(get_job_volume, parameterized_jobs), key=get_job_volume)
    aggregated_jobs_by_volume = groupby(sorted_jobs_by_volume, get_job_volume)
    nomad_pending_jobs_by_volume, nomad_running_jobs_by_volume = _aggregate_nomad_jobs(
        aggregated_jobs_by_volume
    )

    return {
        "nomad_pending_jobs": nomad_pending_jobs,
        "nomad_running_jobs": nomad_running_jobs,
        "nomad_pending_jobs_by_type": nomad_pending_jobs_by_type,
        "nomad_running_jobs_by_type": nomad_running_jobs_by_type,
        "nomad_pending_jobs_by_volume": nomad_pending_jobs_by_volume,
        "nomad_running_jobs_by_volume": nomad_running_jobs_by_volume,
    }


def get_job_details(job):
    """Given a Nomad Job, as returned by the API, returns the type and volume id"""

    # Surveyor jobs don't have ids and RAM, so handle them specially.
    if job["ID"].startswith("SURVEYOR"):
        return "SURVEYOR", False

    # example SALMON_1_2323
    name_match = re.match(r"(?P<type>\w+)_(?P<volume_id>\d+)_\d+$", job["ID"])
    if not name_match:
        return False, False

    return name_match.group("type"), name_match.group("volume_id")


def _aggregate_nomad_jobs(aggregated_jobs):
    """Aggregates the job counts.

    This is accomplished by using the stats that each
    parameterized job has about its children jobs.

    `jobs` should be a response from the Nomad API's jobs endpoint.
    """
    nomad_running_jobs = {}
    nomad_pending_jobs = {}
    for (aggregate_key, group) in aggregated_jobs:
        pending_jobs_count = 0
        running_jobs_count = 0
        for job in group:
            if job["JobSummary"]["Children"]:  # this can be null
                pending_jobs_count += job["JobSummary"]["Children"]["Pending"]
                running_jobs_count += job["JobSummary"]["Children"]["Running"]

        nomad_pending_jobs[aggregate_key] = pending_jobs_count
        nomad_running_jobs[aggregate_key] = running_jobs_count

    return nomad_pending_jobs, nomad_running_jobs


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

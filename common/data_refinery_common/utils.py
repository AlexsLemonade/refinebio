import csv
import os
import requests

from billiard import current_process
from django.core.exceptions import ImproperlyConfigured
from retrying import retry

# Found: http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-metadata.html
METADATA_URL = "http://169.254.169.254/latest/meta-data"
INSTANCE_ID = None


def get_env_variable(var_name: str, default:str=None) -> str:
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

def get_supported_platforms(platforms_csv:str="supported_platforms.csv") -> list:
    """
    Loads our supported platforms file and returns a list of supported
    platform ascession codes.
    CSV must be in the format:
    Species | Platform | Name | Assays | Supported | Processor
    """
    supported_platforms = []
    with open(platforms_csv) as platforms_file:
        reader = csv.reader(platforms_file)
        
        for line in reader:
            
            # Skip the header row
            # Lines are 1 indexed, #BecauseCSV
            if reader.line_num is 1:
                continue

            if line[4] is "Y":
                supported_platforms.append(line[1])

    return supported_platforms
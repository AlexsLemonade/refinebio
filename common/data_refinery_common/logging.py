import boto3
import daiquiri
import logging
import sys
import watchtower

from data_refinery_common.utils import get_worker_id
from data_refinery_common.utils import get_env_variable


# Most of the formatting in this string is for the logging system. All
# that the call to format() does is replace the "{0}" in the string
# with the worker id.
FORMAT_STRING = (
    "%(asctime)s {0} %(name)s %(color)s%(levelname)s%(extras)s"
    ": %(message)s%(color_stop)s"
).format(get_worker_id())


def unconfigure_root_logger():
    root_logger = logging.getLogger(None)

    # Remove all handlers
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)


# We want this function to be run once. Having it top level like this
# will cause it to be run when the module is first imported, then not
# again.
unconfigure_root_logger()


def get_and_configure_logger(name: str) -> logging.Logger:
    # Set level to a environment variable; I think at least
    logger = daiquiri.getLogger(name)
    logger.setLevel(logging.INFO)

    # This is the local handler
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(daiquiri.formatter.ColorExtrasFormatter(
        fmt=FORMAT_STRING, keywords=[]))
    logger.logger.addHandler(handler)

    # This is for CloudWatching.
    # Watchtower honestly isn't an awesome package and 
    # we should keep an eye on it and possibly reevaluate if
    # it becomes problematic.
    # Related: https://github.com/kislyuk/watchtower/issues/34
    # Related: https://github.com/kislyuk/watchtower/issues/56
    in_cloud = get_env_variable("RUNNING_IN_CLOUD", False)
    # Env vars are only ever strings.
    if in_cloud == "True":
        # <service>-<user>-<stage>
        user = get_env_variable("USER", "user")
        stage = get_env_variable("STAGE", "dev")
        service = get_env_variable("SERVICE", "common")

        # AWS
        aws_access_key_id = get_env_variable("aws_access_key_id".upper(), "AK123")
        aws_secret_access_key = get_env_variable("aws_secret_access_key".upper(), "SK123")
        region = get_env_variable("REGION".upper(), "us-east-1")

        log_group = "data-refinery-log-group-{user}-{stage}".format(user=user, stage=stage)
        stream_name = "log-stream-{service}-django-{user}-{stage}".format(service=service, user=user, stage=stage)
        session = boto3.session.Session(aws_access_key_id=aws_access_key_id, 
            aws_secret_access_key=aws_secret_access_key, 
            region_name=region)

        try:
            tower = watchtower.CloudWatchLogHandler(
                log_group=log_group,
                stream_name=stream_name,
                boto3_session=session,
                create_log_group=False
                )
            logger.logger.addHandler(tower)
        except Exception as e:
            # Oh, the irony!
            logger.error("Error while creating CloudwatchLogHandler.")

    return logger

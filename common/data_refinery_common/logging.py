import boto3
import daiquiri
import logging
import sys

from data_refinery_common.utils import get_worker_id
from data_refinery_common.utils import get_env_variable_gracefully


# Most of the formatting in this string is for the logging system. All
# that the call to format() does is replace the "{0}" in the string
# with the worker id.
FORMAT_STRING = (
    "%(asctime)s {0} %(name)s %(color)s%(levelname)s%(extras)s"
    ": %(message)s%(color_stop)s"
).format(get_worker_id())


def unconfigure_root_logger():
    """Prevents the root logger from duplicating our messages.

    The root handler comes preconfigured with a handler. This causes
    all our logs to be logged twice, once with our cool handler and
    one that lacks all context. This function removes that stupid
    extra handler.
    """
    root_logger = logging.getLogger(None)
    # Remove all handlers
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)


def get_and_configure_logger(name: str) -> logging.Logger:
    unconfigure_root_logger()

    # Set level to a environment variable; I think at least
    logger = daiquiri.getLogger(name)
    logger.setLevel(logging.INFO)

    # This is the local handler
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(daiquiri.formatter.ColorExtrasFormatter(
        fmt=FORMAT_STRING, keywords=[]))
    logger.logger.addHandler(handler)

    # This is the Sentry handler
    if "data_refinery_api" in name:
        raven_dsn = get_env_variable_gracefully("RAVEN_DSN_API", False)
    else:
        raven_dsn = get_env_variable_gracefully("RAVEN_DSN", False)
    if raven_dsn:
        from raven.contrib.django.handlers import SentryHandler

        handler = SentryHandler()
        handler.setFormatter(daiquiri.formatter.ColorExtrasFormatter(
            fmt=FORMAT_STRING, keywords=[]))
        handler.setLevel(logging.WARNING)
        logger.logger.addHandler(handler)

    return logger

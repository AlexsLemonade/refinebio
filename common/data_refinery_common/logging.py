from typing import List
import daiquiri
import logging
import sys
from data_refinery_common.utils import get_worker_id


base_format_string = "%(asctime)s {0} %(name)s %(color)s%(levelname)s"
parameter_format_string = " [{0}: %({1})s]"


def unconfigure_root_logger():
    root_logger = logging.getLogger(None)

    # Remove all handlers
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)


def get_logger_with_parameters(name: str, parameters: List[str]) -> logging.Logger:
    unconfigure_root_logger()
    # Set level to a environment variable
    logger = daiquiri.getLogger(name, level=logging.INFO)
    logger.propagate = False
    format_string = base_format_string.format(get_worker_id())
    for parameter in parameters:
        format_string = format_string + parameter_format_string.format(parameter, parameter)

    format_string = format_string + ": %(message)s%(color_stop)s"

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(daiquiri.formatter.ColorFormatter(fmt=format_string))
    logger.logger.addHandler(handler)
    return logger

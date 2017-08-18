from typing import List
import daiquiri
import logging
import sys
import weakref
from data_refinery_common.utils import get_worker_id


# base_format_string = "%(asctime)s {0} %(name)s %(color)s%(levelname)s"
format_string = (
    "%(asctime)s {0} %(name)s %(color)s%(levelname)s%(context)s"
    ": %(message)s%(color_stop)s"
).format(get_worker_id())
parameter_template = " [{0}: {1}]"


class ArbitraryContextAdapter(logging.LoggerAdapter):
    """Logger adapter to add context string to log record's extra data

    Keywords passed to the log call are formatted into a context
    string, which is then added to the "extra" dictionary passed to
    the underlying logger so they are emitted with the log message and
    available to the format string.

    Example:
      A formatter with the format string of "%(context)s %(messages)"
      called like logger.info("A message.", test1="a", test2="b")
      would log the message: " [test1: a] [test2: b] A message."

    Special keywords:

    extra
      An existing dictionary of extra values to be passed to the
      logger. If present, the dictionary is copied and extended.
    """

    def process(self, msg, kwargs):
        extra = self.extra.copy()
        if "extra" in kwargs:
            extra.update(kwargs.pop("extra"))

        context = ""
        for name in list(kwargs.keys()):
            if name == "exc_info":
                continue
            context += parameter_template.format(name, kwargs.pop(name))

        extra["context"] = context
        extra["_daiquiri_extra"] = extra
        kwargs["extra"] = extra
        return msg, kwargs


_LOGGERS = weakref.WeakValueDictionary()


def getLogger(name=None, **kwargs):
    """Build a logger with the given name.
    :param name: The name for the logger. This is usually the module
                 name, ``__name__``.
    :type name: string
    """
    adapter = _LOGGERS.get(name)
    if not adapter:
        # The following note was found in the original source found here:
        # https://github.com/jd/daiquiri/blob/master/daiquiri/__init__.py
        # NOTE(jd) Keep using the `adapter' variable here because so it's not
        # collected by Python since _LOGGERS contains only a weakref
        adapter = ArbitraryContextAdapter(logging.getLogger(name), kwargs)
        _LOGGERS[name] = adapter
    return adapter


def unconfigure_root_logger():
    root_logger = logging.getLogger(None)

    # Remove all handlers
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)


def get_and_configure_logger(name: str) -> logging.Logger:
    unconfigure_root_logger()
    # Set level to a environment variable; I think at least
    logger = getLogger(name, level=logging.INFO)

    # I think using both this and unconfigure_root_logger might be redundant
    # logger.propagate = False

    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(daiquiri.formatter.ColorFormatter(fmt=format_string))
    logger.logger.addHandler(handler)
    return logger

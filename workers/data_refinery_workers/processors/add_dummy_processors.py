"""Create dummy processor records that are defined in ProcessorEnum.
This module will be invoked by all test modules.
"""

from data_refinery_common.models import Processor
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import _names

for label in _names.ProcessorEnum:
    Processor.objects.get_or_create(name=label.value, version=__version__)

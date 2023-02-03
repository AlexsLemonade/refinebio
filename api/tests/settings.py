"""
Django settings for data_refinery_api project tests.
"""

from data_refinery_api.settings import *

ELASTICSEARCH_INDEX_NAMES = {
    "data_refinery_common.models.documents": "experiments_test",
}
ELASTICSEARCH_DSL_AUTOSYNC = True

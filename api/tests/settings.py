"""
Django settings for data_refinery_api project tests.
"""

from data_refinery_api.settings import *

OPENSEARCH_INDEX_NAMES = {
    "data_refinery_common.models.documents": "experiments_test",
}
OPENSEARCH_DSL_AUTOSYNC = True

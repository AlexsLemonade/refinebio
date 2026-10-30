"""
Django settings for data_refinery_api project tests.
"""

from data_refinery_api.settings import *

OPENSEARCH_INDEX_NAMES = {
    "data_refinery_common.models.documents": "experiments_test",
}
OPENSEARCH_DSL_AUTOSYNC = True
# Force an immediate refresh on auto-sync writes so tests see newly-saved
# docs on the next query. Prod leaves this False (the package default) to
# avoid the per-write refresh cost.
OPENSEARCH_DSL_AUTO_REFRESH = True

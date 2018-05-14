from __future__ import absolute_import, unicode_literals

import os
import string
import warnings
from django.utils import timezone
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OriginalFile, ComputationalResult, ComputedFile, SampleResultAssociation
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.utils import get_env_variable
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """
    Also adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """

    return job_context

def _smash(job_context: Dict) -> Dict:
    """
    Smash all of the samples together!

    Steps:
        Combine common genes (Pandas merge? ^)
        Transpose such that genes are columns (features)
        Scale features with sci-kit learn
        Transpose again such that samples are columns and genes are rows
    """

    return job_context

def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    return job_context

def smash(job_id: int) -> None:
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _determine_brainarray_package,
                        _run_scan_upc,
                        _create_result_objects,
                        # utils.upload_processed_files,
                        # utils.cleanup_raw_files,
                        utils.end_job])

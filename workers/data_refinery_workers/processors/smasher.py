from __future__ import absolute_import, unicode_literals

import os
import string
import warnings
from django.utils import timezone
from typing import Dict

import pandas as pd
from sklearn import preprocessing

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OriginalFile, ComputationalResult, ComputedFile, SampleResultAssociation
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.utils import get_env_variable
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """
    Fetches and prepares the files to smash.
    """

    all_sample_files = []
    for sample in job_context["samples"]:
        all_sample_files = all_sample_files + list(sample.get_result_files())
    all_sample_files = list(set(all_sample_files))
    job_context['input_files'] = all_sample_files

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

    # Merge all the frames into one
    all_frames = []
    for computed_file in job_context['input_files']:
        data = pd.DataFrame.from_csv(computed_file.absolute_file_path, sep='\t', header=0)
        all_frames.append(data)
    merged = all_frames[0]
    i = 1
    while i < len(all_frames):
        merged = merged.merge(all_frames[i], left_index=True, right_index=True)
        i = i + 1

    # Transpose before scaling
    transposed = merged.transpose()
    
    # Scale
    # XXX/TODO: Is MinMaxScaler or Standard or Robust the right scaler here?
    scaler = preprocessing.MinMaxScaler(copy=True)
    scaler.fit(transposed)
    scaled = pd.DataFrame(scaler.transform(transposed), index=transposed.index, columns=transposed.columns)

    # Untranspose
    untransposed = scaled.transpose()

    # Write to temp file
    smash_path = "/home/user/data_store/smashed/"
    outfile = smash_path + str(job_context["dataset"].pk) + ".csv"
    os.makedirs(smash_path, exist_ok=True)
    untransposed.to_csv(outfile, sep='\t', encoding='utf-8')
    job_context["output_file"] = outfile

    return job_context

def _upload_and_notify(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """
    return job_context

def _update_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a run is complete """

    dataset = job_context["dataset"]
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.save()

    return job_context

def smash(job_id: int) -> None:

    return utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _prepare_files,
                        _smash,
                        _upload_and_notify,
                        _update_result_objects,
                        utils.end_job])

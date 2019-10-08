# -*- coding: utf-8 -*-

import boto3
import csv
import os
import rpy2
import rpy2.robjects as ro
import shutil
import simplejson as json
import string
import warnings
import requests
import psutil
import multiprocessing
import logging
import time

from botocore.exceptions import ClientError
from datetime import datetime, timedelta
from django.conf import settings
from django.utils import timezone
from pathlib import Path
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr
from sklearn import preprocessing
from typing import Dict, List
import numpy as np
import pandas as pd

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    OriginalFile,
    Pipeline,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable, calculate_file_size, calculate_sha1
from data_refinery_workers.processors import utils
from urllib.parse import quote


RESULTS_BUCKET = get_env_variable("S3_RESULTS_BUCKET_NAME", "refinebio-results-bucket")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
BODY_HTML = Path('data_refinery_workers/processors/smasher_email.min.html').read_text().replace('\n', '')
BODY_ERROR_HTML = Path('data_refinery_workers/processors/smasher_email_error.min.html').read_text().replace('\n', '')
BYTES_IN_GB = 1024 * 1024 * 1024
logger = get_and_configure_logger(__name__)
### DEBUG ###
logger.setLevel(logging.getLevelName('DEBUG'))



def log_failure(job_context: Dict, failure_reason: str) -> Dict:
    if not job_context["job"].failure_reason:
        job_context["job"].failure_reason = failure_reason

    job_context['success'] = False
    logger.warn(job_context["job"].failure_reason, job_id=job_context['job'].id)

    return job_context


def log_state(message, job, start_time=False):
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process(os.getpid())
        ram_in_GB = process.memory_info().rss / BYTES_IN_GB
        logger.debug(message,
                     total_cpu=psutil.cpu_percent(),
                     process_ram=ram_in_GB,
                     job_id=job.id)

        if start_time:
            logger.debug('Duration: %s' % (time.time() - start_time), job_id=job.id)
        else:
            return time.time()


def prepare_files(job_context: Dict) -> Dict:
    """
    Fetches and prepares the files to smash.
    """
    start_prepare_files = log_state("start prepare files", job_context["job"])
    found_files = False
    job_context['input_files'] = {}
    # `key` can either be the species name or experiment accession.
    for key, samples in job_context["samples"].items():
        smashable_files = []
        seen_files = set()
        for sample in samples:
            smashable_file = sample.get_most_recent_smashable_result_file()
            if smashable_file is not None and smashable_file not in seen_files:
                smashable_files = smashable_files + [(smashable_file, sample)]
                seen_files.add(smashable_file)
                found_files = True

        job_context['input_files'][key] = smashable_files

    if not found_files:
        error_message = "Couldn't get any files to smash for Smash job!!"
        logger.error(error_message,
                     dataset_id=job_context['dataset'].id,
                     num_samples=len(job_context["samples"]))

        # Delay failing this pipeline until the failure notify has been sent
        job_context['dataset'].failure_reason = error_message
        job_context['dataset'].success = False
        job_context['dataset'].save()
        job_context['job'].success = False
        job_context["job"].failure_reason = "Couldn't get any files to smash for Smash job - empty all_sample_files"
        return job_context

    job_context["work_dir"] = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk) + "/"
    # Ensure we have a fresh smash directory
    shutil.rmtree(job_context["work_dir"], ignore_errors=True)
    os.makedirs(job_context["work_dir"])

    job_context["output_dir"] = job_context["work_dir"] + "output/"
    os.makedirs(job_context["output_dir"])
    log_state("end prepare files", job_context["job"], start_prepare_files)
    return job_context


def _load_and_sanitize_file(computed_file_path):
    """ Read and sanitize a computed file """

    data = pd.read_csv(computed_file_path, sep='\t', header=0, index_col=0, error_bad_lines=False)

    # Strip any funky whitespace
    data.columns = data.columns.str.strip()
    data = data.dropna(axis='columns', how='all')

    # Make sure the index type is correct
    data.index = data.index.map(str)

    # Ensure that we don't have any dangling Brainarray-generated probe symbols.
    # BA likes to leave '_at', signifying probe identifiers,
    # on their converted, non-probe identifiers. It makes no sense.
    # So, we chop them off and don't worry about it.
    data.index = data.index.str.replace('_at', '')

    # Remove any lingering Affymetrix control probes ("AFFX-")
    data = data[~data.index.str.contains('AFFX-')]

    # If there are any _versioned_ gene identifiers, remove that
    # version information. We're using the latest brainarray for everything anyway.
    # Jackie says this is okay.
    # She also says that in the future, we may only want to do this
    # for cross-technology smashes.

    # This regex needs to be able to handle EGIDs in the form:
    #       ENSGXXXXYYYZZZZ.6
    # and
    #       fgenesh2_kg.7__3016__AT5G35080.1 (via http://plants.ensembl.org/Arabidopsis_lyrata/Gene/Summary?g=fgenesh2_kg.7__3016__AT5G35080.1;r=7:17949732-17952000;t=fgenesh2_kg.7__3016__AT5G35080.1;db=core)
    data.index = data.index.str.replace(r"(\.[^.]*)$", '')

    # Squish duplicated rows together.
    # XXX/TODO: Is mean the appropriate method here?
    #           We can make this an option in future.
    # Discussion here: https://github.com/AlexsLemonade/refinebio/issues/186#issuecomment-395516419
    data = data.groupby(data.index, sort=False).mean()

    return data


def process_frame(inputs) -> Dict:
    (job_context, computed_file, sample, index) = inputs
    log_state('processing frame {}'.format(index), job_context["job"])
    frame = {
        "unsmashable": False,
        "unsmashable_file": None,
        "technology": None,
        "dataframe": None
    }

    def unsmashable(path):
      frame['unsmashable'] = True
      frame['unsmashable_file'] = path
      return frame

    def smashable(data):
      frame['dataframe'] = data
      return frame

    try:
        # Download the file to a job-specific location so it
        # won't disappear while we're using it.
        computed_file_path = job_context["work_dir"] + computed_file.filename
        computed_file_path = computed_file.get_synced_file_path(path=computed_file_path)

        # Bail appropriately if this isn't a real file.
        if not computed_file_path or not os.path.exists(computed_file_path):
            logger.warning("Smasher received non-existent file path.",
                computed_file_path=computed_file_path,
                computed_file=computed_file,
                dataset_id=job_context['dataset'].id,
            )
            return unsmashable(computed_file_path)

        data = _load_and_sanitize_file(computed_file_path)

        if len(data.columns) > 2:
            # Most of the time, >1 is actually bad, but we also need to support
            # two-channel samples. I think ultimately those should be given some kind of
            # special consideration.
            logger.info("Found a frame with more than 2 columns - this shouldn't happen!",
                computed_file_path=computed_file_path,
                computed_file_id=computed_file.id
            )
            return unsmashable(computed_file_path)

        # via https://github.com/AlexsLemonade/refinebio/issues/330:
        #   aggregating by experiment -> return untransformed output from tximport
        #   aggregating by species -> log2(x + 1) tximport output
        if job_context['dataset'].aggregate_by == 'SPECIES' \
        and computed_file_path.endswith("lengthScaledTPM.tsv"):
            data = data + 1
            data = np.log2(data)

        # Detect if this data hasn't been log2 scaled yet.
        # Ideally done in the NO-OPPER, but sanity check here.
        if (not computed_file_path.endswith("lengthScaledTPM.tsv")) and (data.max() > 100).any():
            logger.info("Detected non-log2 microarray data.", file=computed_file)
            data = np.log2(data)

        # Explicitly title this dataframe
        try:
            # Unfortuantely, we can't use this as `title` can cause a collision
            # data.columns = [computed_file.samples.all()[0].title]
            # So we use this, which also helps us support the case of missing SampleComputedFileAssociation
            data.columns = [sample.accession_code]
        except ValueError as e:
            # This sample might have multiple channels, or something else.
            # Don't mess with it.
            logger.warn("Smasher found multi-channel column (probably) - skipping!",
                        exc_info=1,
                        computed_file_path=computed_file_path,
            )
            return unsmashable(computed_file.filename)
        except Exception as e:
            # Okay, somebody probably forgot to create a SampleComputedFileAssociation
            # Don't mess with it.
            logger.warn("Smasher found very bad column title - skipping!",
                        exc_info=1,
                        computed_file_path=computed_file_path
            )
            return unsmashable(computed_file.filename)

        is_rnaseq = computed_file_path.endswith("lengthScaledTPM.tsv")
        frame['technology'] = 'rnaseq' if is_rnaseq else 'microarray'

    except Exception as e:
        logger.exception("Unable to smash file",
            file=computed_file_path,
            dataset_id=job_context['dataset'].id,
        )
        return unsmashable(computed_file_path)
    finally:
        # Delete before archiving the work dir
        if computed_file_path and os.path.exists(computed_file_path):
            os.remove(computed_file_path)

    return smashable(data)


def process_frames_for_key(key: str, input_files: List[ComputedFile], job_context: Dict) -> Dict:
    job_context['original_merged'] = pd.DataFrame()

    start_frames = log_state("building frames for species or experiment {}".format(key), job_context["job"])
    # Merge all the frames into one
    #cpus = max(1, psutil.cpu_count()/2)
    #with multiprocessing.Pool(int(cpus)) as pool:
    #    processed_frames = pool.map(
    mapped_frames = map(
        process_frame,
        [(
            job_context,
            computed_file,
            sample,
            i
        ) for i, (computed_file, sample) in enumerate(input_files)
    ])
    processed_frames = list(mapped_frames)

    # Build up a list of microarray frames and a list of
    # rnaseq frames and then combine them so they're sorted
    # out.
    job_context['microarray_frames'] = []
    job_context['rnaseq_frames'] = []
    for frame in processed_frames:
        if frame['technology'] is 'microarray':
            job_context['microarray_frames'].append(frame['dataframe'])
        elif frame['technology'] is 'rnaseq':
            job_context['rnaseq_frames'].append(frame['dataframe'])

        if frame['unsmashable']:
            job_context['unsmashable_files'].append(frame['unsmashable_file'])

    job_context['num_samples'] = job_context['num_samples'] + len(job_context['microarray_frames'])
    job_context['num_samples'] = job_context['num_samples'] + len(job_context['rnaseq_frames'])

    log_state("set frames for key {}".format(key), job_context["job"], start_frames)

    return job_context


def compile_metadata(job_context: Dict) -> Dict:
    """Compiles metadata about the job.

    Returns a new dict containing the metadata, not the job_context.
    """
    if not job_context['dataset'].quant_sf_only:
        metadata['aggregate_by'] = job_context["dataset"].aggregate_by
        metadata['scale_by'] = job_context["dataset"].scale_by
        # https://github.com/AlexsLemonade/refinebio/pull/421#discussion_r203799646
        metadata['non_aggregated_files'] = unsmashable_files
        metadata['ks_statistic'] = job_context.get("ks_statistic", None)
        metadata['ks_pvalue'] = job_context.get("ks_pvalue", None)
        metadata['ks_warning'] = job_context.get("ks_warning", None)
        metadata['quantile_normalized'] = job_context['dataset'].quantile_normalize

    samples = {}
    for sample in job_context["dataset"].get_samples():
        samples[sample.accession_code] = sample.to_metadata_dict()

    metadata['samples'] = samples

    experiments = {}
    for experiment in job_context["dataset"].get_experiments():
        exp_dict = experiment.to_metadata_dict()
        exp_dict['sample_accession_codes'] = [v for v in experiment.samples.all().values_list('accession_code', flat=True)]
        experiments[experiment.accession_code] = exp_dict

    metadata['experiments'] = experiments

    return metadata

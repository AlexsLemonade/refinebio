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
import logging
import time
from multiprocessing import Pool

from botocore.exceptions import ClientError
from datetime import timedelta
from django.conf import settings
from django.utils import timezone
from pathlib import Path
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr
from sklearn import preprocessing
from typing import Dict, List, Tuple
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
from data_refinery_workers.processors import utils, smashing_utils
from urllib.parse import quote


RESULTS_BUCKET = get_env_variable("S3_RESULTS_BUCKET_NAME", "refinebio-results-bucket")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
BODY_HTML = Path('data_refinery_workers/processors/smasher_email.min.html').read_text().replace('\n', '')
BODY_ERROR_HTML = Path('data_refinery_workers/processors/smasher_email_error.min.html').read_text().replace('\n', '')
BYTES_IN_GB = 1024 * 1024 * 1024
logger = get_and_configure_logger(__name__)
### DEBUG ###
logger.setLevel(logging.getLevelName('DEBUG'))

PROCESS_POOL_SIZE = max(1, int(psutil.cpu_count()/2 - 1))

SCALERS = {
    'MINMAX': preprocessing.MinMaxScaler,
    'STANDARD': preprocessing.StandardScaler,
    'ROBUST': preprocessing.RobustScaler,
}


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


def _prepare_files(job_context: Dict) -> Dict:
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

def downlad_computed_file(download_tuple: Tuple[ComputedFile, str]):
    """ this function downloads the latest computed file. Receives a tuple with
    the computed file and the path where it needs to be downloaded
    This is used to parallelize downloading quantsf files. """
    (latest_computed_file, output_file_path) = download_tuple
    try:
        latest_computed_file.get_synced_file_path(path=output_file_path)
    except:
        # Let's not fail if there's an error syncing one of the quant.sf files
        logger.exception('Failed to sync computed file', computed_file_id=latest_computed_file.pk)

def sync_quant_files(output_path, files_sample_tuple):
    """ Takes a list of ComputedFiles and copies the ones that are quant files to the provided directory.
        Returns the total number of samples that were included """
    num_samples = 0
    samples = [sample for (_, sample) in files_sample_tuple]
    page_size = 100
    # split the samples in groups and download each one individually
    pool = Pool(processes=PROCESS_POOL_SIZE)

    # for each sample we need it's latest quant.sf file we don't want to query the db
    # for all of them, so we do it in groups of 100, and then download all of the computed_files
    # in parallel
    for sample_page in (samples[i*page_size:i+page_size] for i in range(0, len(samples), page_size)):
        sample_and_computed_files = []
        for sample in sample_page:
            latest_computed_file = sample.get_most_recent_quant_sf_file()
            if not latest_computed_file:
                continue
            output_file_path = output_path + sample.accession_code + "_quant.sf"
            sample_and_computed_files.append((latest_computed_file, output_file_path))

        # download this set of files, this will take a few seconds that should also help the db recover
        pool.map(downlad_computed_file, sample_and_computed_files)
        num_samples += len(sample_and_computed_files)

    pool.close()
    pool.join()

    return num_samples

def _inner_join(job_context: Dict) -> pd.DataFrame:
    """Performs an inner join across the all_frames key of job_context.

    Returns a new dict containing the metadata, not the job_context.
    """
    # Merge all of the frames we've gathered into a single big frame, skipping duplicates.
    # TODO: If the very first frame is the wrong platform, are we boned?
    merged = job_context['all_frames'][0]
    i = 1

    old_len_merged = len(merged)
    merged_backup = merged

    while i < len(job_context['all_frames']):
        frame = job_context['all_frames'][i]
        i = i + 1

        if i % 1000 == 0:
            logger.info("Smashing keyframe",
                        i=i,
                        job_id=job_context['job'].id)

        # I'm not sure where these are sneaking in from, but we don't want them.
        # Related: https://github.com/AlexsLemonade/refinebio/issues/390
        breaker = False
        for column in frame.columns:
            if column in merged.columns:
                breaker = True

        if breaker:
            logger.warning("Column repeated for smash job!",
                           input_files=str(input_files),
                           dataset_id=job_context["dataset"].id,
                           job_id=job_context["job"].id,
                           column=column)
            continue

        # This is the inner join, the main "Smash" operation
        merged = merged.merge(frame, how='inner', left_index=True, right_index=True)

        new_len_merged = len(merged)
        if new_len_merged < old_len_merged:
            logger.warning("Dropped rows while smashing!",
                dataset_id=job_context["dataset"].id,
                old_len_merged=old_len_merged,
                new_len_merged=new_len_merged
            )
        if new_len_merged == 0:
            logger.warning("Skipping a bad merge frame!",
                           dataset_id=job_context["dataset"].id,
                           job_id=job_context["job"].id,
                           old_len_merged=old_len_merged,
                           new_len_merged=new_len_merged,
                           bad_frame_number=i,)
            merged = merged_backup
            new_len_merged = len(merged)
            try:
                job_context['unsmashable_files'].append(frame.columns[0])
            except Exception:
                # Something is really, really wrong with this frame.
                pass

        old_len_merged = len(merged)
        merged_backup = merged

    return merged


def _smash_key(job_context: Dict, key: str, input_files: List[ComputedFile]) -> Dict:
    """Smash all of the input files together for a given key.

    Steps:
        Combine common genes (pandas merge)
        Transpose such that genes are columns (features)
        Scale features with sci-kit learn
        Transpose again such that samples are columns and genes are rows
    """
    start_smash = log_state("start _smash_key for {}".format(key), job_context["job"])

    # Check if we need to copy the quant.sf files
    if job_context['dataset'].quant_sf_only:
        outfile_dir = job_context["output_dir"] + key + "/"
        os.makedirs(outfile_dir, exist_ok=True)
        job_context['num_samples'] += sync_quant_files(outfile_dir, input_files)
        # we ONLY want to give quant sf files to the user if that's what they requested
        return job_context

    job_context = smashing_utils.process_frames_for_key(key,
                                                        input_files,
                                                        job_context,
                                                        merge_strategy='inner')

    # Combine the two technologies into a single list of dataframes.
    ## Extend one list rather than adding the two together so we don't
    ## the memory both are using.
    ## Also free up the the memory the microarray-only list was using with pop.
    job_context['rnaseq_frames'].extend(job_context.pop('microarray_frames'))
    ## Change the key of the now-extended list
    job_context['all_frames'] = job_context.pop('rnaseq_frames')

    if len(job_context['all_frames']) < 1:
        logger.error("Was told to smash a key with no frames!",
                     job_id=job_context['job'].id,
                     key=key)
        # TODO: is this the proper way to handle this? I can see us
        # not wanting to fail an entire dataset because one experiment
        # had a problem, but I also think it could be problematic to
        # just skip an experiment and pretend nothing went wrong.
        return job_context

    merged = _inner_join(job_context)

    job_context['original_merged'] = merged
    log_state("end build all frames", job_context["job"], start_smash)
    start_qn = log_state("start qn", job_context["job"], start_smash)

    # Quantile Normalization
    if job_context['dataset'].quantile_normalize:
        try:
            job_context['merged_no_qn'] = merged
            job_context['organism'] = job_context['dataset'].get_samples().first().organism
            job_context = smashing_utils.quantile_normalize(job_context)
            merged = job_context.get('merged_qn', None)

            # We probably don't have an QN target or there is another error,
            # so let's fail gracefully.
            assert merged is not None, "Problem occured during quantile normalization: No merged_qn"
        except Exception as e:
            logger.exception("Problem occured during quantile normalization",
                dataset_id=job_context['dataset'].id,
                processor_job_id=job_context["job"].id,
            )
            job_context['dataset'].success = False

            if not job_context['job'].failure_reason:
                job_context['job'].failure_reason = "Failure reason: " + str(e)
                job_context['dataset'].failure_reason = "Failure reason: " + str(e)

            job_context['dataset'].save()
            # Delay failing this pipeline until the failure notify has been sent
            job_context['job'].success = False
            job_context['failure_reason'] = str(e)
            return job_context

    # End QN
    log_state("end qn", job_context["job"], start_qn)
    # Transpose before scaling
    # Do this even if we don't want to scale in case transpose
    # modifies the data in any way. (Which it shouldn't but
    # we're paranoid.)
    # TODO: stop the paranoia because Josh has alleviated it.
    transposed = merged.transpose()
    start_scaler = log_state("starting scaler", job_context["job"])
    # Scaler
    if job_context['dataset'].scale_by != "NONE":
        scale_funtion = SCALERS[job_context['dataset'].scale_by]
        scaler = scale_funtion(copy=True)
        scaler.fit(transposed)
        scaled = pd.DataFrame(  scaler.transform(transposed),
                                index=transposed.index,
                                columns=transposed.columns
                            )
        # Untranspose
        untransposed = scaled.transpose()
    else:
        # Wheeeeeeeeeee
        untransposed = transposed.transpose()
    log_state("end scaler", job_context["job"], start_scaler)

    # This is just for quality assurance in tests.
    job_context['final_frame'] = untransposed

    # Write to temp file with dataset UUID in filename.
    subdir = ''
    if job_context['dataset'].aggregate_by in ["SPECIES", "EXPERIMENT"]:
        subdir = key
    elif job_context['dataset'].aggregate_by == "ALL":
        subdir = "ALL"

    # Normalize the Header format
    untransposed.index.rename('Gene', inplace=True)

    outfile_dir = job_context["output_dir"] + key + "/"
    os.makedirs(outfile_dir, exist_ok=True)
    outfile = outfile_dir + key + ".tsv"
    job_context['smash_outfile'] = outfile
    untransposed.to_csv(outfile, sep='\t', encoding='utf-8')

    log_state("end _smash_key for {}".format(key), job_context["job"], start_smash)

    return job_context


def _smash_all(job_context: Dict) -> Dict:
    """Perform smashing on all species/experiments in the dataset.
    """
    start_smash = log_state("start smash", job_context["job"])
    # We have already failed - return now so we can send our fail email.
    if job_context['job'].success is False:
        return job_context

    try:
        job_context['unsmashable_files'] = []
        job_context['num_samples'] = 0

        # Smash all of the sample sets
        logger.debug("About to smash!",
                     dataset_count=len(job_context['dataset'].data),
                     job_id=job_context['job'].id)

        # Once again, `key` is either a species name or an experiment accession
        for key, input_files in job_context['input_files'].items():
            job_context = _smash_key(job_context, key, input_files)

        smashing_utils.write_non_data_files(job_context)

        # Finally, compress all files into a zip
        final_zip_base = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk)
        shutil.make_archive(final_zip_base, 'zip', job_context["output_dir"])
        job_context["output_file"] = final_zip_base + ".zip"
    except Exception as e:
        logger.exception("Could not smash dataset.",
                        dataset_id=job_context['dataset'].id,
                        processor_job_id=job_context['job_id'],
                        num_input_files=len(job_context['input_files']))
        job_context['dataset'].success = False
        job_context['job'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].save()
        # Delay failing this pipeline until the failure notify has been sent
        job_context['job'].success = False
        job_context['failure_reason'] = str(e)
        return job_context

    job_context['dataset'].success = True
    job_context['dataset'].save()

    logger.debug("Created smash output!",
        archive_location=job_context["output_file"])

    log_state("end smash", job_context["job"], start_smash);
    return job_context


def _upload(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """

    # There has been a failure already, don't try to upload anything.
    if not job_context.get("output_file", None):
        logger.error("Was told to upload a smash result without an output_file.",
                     job_id=job_context['job'].id)
        return job_context

    try:
        if job_context.get("upload", True) and settings.RUNNING_IN_CLOUD:
            s3_client = boto3.client('s3')

            # Note that file expiry is handled by the S3 object lifecycle,
            # managed by terraform.
            s3_client.upload_file(
                    job_context["output_file"],
                    RESULTS_BUCKET,
                    job_context["output_file"].split('/')[-1],
                    ExtraArgs={'ACL':'public-read'}
                )
            result_url = ("https://s3.amazonaws.com/" + RESULTS_BUCKET + "/" +
                          job_context["output_file"].split('/')[-1])

            job_context["result_url"] = result_url

            logger.debug("Result uploaded!",
                    result_url=job_context["result_url"]
                )

            job_context["dataset"].s3_bucket = RESULTS_BUCKET
            job_context["dataset"].s3_key = job_context["output_file"].split('/')[-1]
            job_context["dataset"].size_in_bytes = calculate_file_size(job_context["output_file"])
            job_context["dataset"].sha1 = calculate_sha1(job_context["output_file"])

            job_context["dataset"].save()

            # File is uploaded, we can delete the local.
            try:
                os.remove(job_context["output_file"])
            except OSError:
                pass

    except Exception as e:
        logger.exception("Failed to upload smash result file.", file=job_context["output_file"])
        job_context['job'].success = False
        job_context['job'].failure_reason = "Failure reason: " + str(e)
        # Delay failing this pipeline until the failure notify has been sent
        # job_context['success'] = False

    return job_context

def _notify(job_context: Dict) -> Dict:
    """ Use AWS SES to notify a user of a smash result.. """

    ##
    # SES
    ##
    if job_context.get("upload", True) and settings.RUNNING_IN_CLOUD:
        # Link to the dataset page, where the user can re-try the download job
        dataset_url = 'https://www.refine.bio/dataset/' + str(job_context['dataset'].id)

        # Send a notification to slack when a dataset fails to be processed
        if job_context['job'].success is False:
            try:
                requests.post(
                    "https://hooks.slack.com/services/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",
                    json={
                        'fallback': 'Dataset failed processing.',
                        'title': 'Dataset failed processing',
                        'title_link': dataset_url,
                        "attachments":[
                            {
                                "color": "warning",
                                "text": job_context['job'].failure_reason,
                                'author_name': job_context["dataset"].email_address,
                                'fields': [
                                    {
                                        'title': 'Dataset id',
                                        'value': str(job_context['dataset'].id)
                                    }
                                ]
                            }
                        ]
                    },
                    headers={'Content-Type': 'application/json'},
                    timeout=10
                )
            except Exception as e:
                logger.warn(e) # It doens't really matter if this didn't work
                pass

        # Don't send an email if we don't have address.
        if job_context["dataset"].email_address:
            SENDER = "Refine.bio Mail Robot <noreply@refine.bio>"
            RECIPIENT = job_context["dataset"].email_address
            AWS_REGION = "us-east-1"
            CHARSET = "UTF-8"


            if job_context['job'].success is False:
                SUBJECT = "There was a problem processing your refine.bio dataset :("
                BODY_TEXT = "We tried but were unable to process your requested dataset. Error was: \n\n" + str(job_context['job'].failure_reason) + "\nDataset ID: " + str(job_context['dataset'].id) + "\n We have been notified and are looking into the problem. \n\nSorry!"

                ERROR_EMAIL_TITLE = quote('I can\'t download my dataset')
                ERROR_EMAIL_BODY = quote("""
                [What browser are you using?]
                [Add details of the issue you are facing]

                ---
                """ + str(job_context['dataset'].id))

                FORMATTED_HTML = BODY_ERROR_HTML.replace('REPLACE_DATASET_URL', dataset_url)\
                                                .replace('REPLACE_ERROR_TEXT', job_context['job'].failure_reason)\
                                                .replace('REPLACE_NEW_ISSUE', 'https://github.com/AlexsLemonade/refinebio/issues/new?title={0}&body={1}&labels=bug'.format(ERROR_EMAIL_TITLE, ERROR_EMAIL_BODY))\
                                                .replace('REPLACE_MAILTO', 'mailto:ccdl@alexslemonade.org?subject={0}&body={1}'.format(ERROR_EMAIL_TITLE, ERROR_EMAIL_BODY))
                job_context['success'] = False
            else:
                SUBJECT = "Your refine.bio Dataset is Ready!"
                BODY_TEXT = "Hot off the presses:\n\n" + dataset_url + "\n\nLove!,\nThe refine.bio Team"
                FORMATTED_HTML = BODY_HTML.replace('REPLACE_DOWNLOAD_URL', dataset_url)\
                                          .replace('REPLACE_DATASET_URL', dataset_url)

            # Try to send the email.
            try:

                # Create a new SES resource and specify a region.
                client = boto3.client('ses', region_name=AWS_REGION)

                #Provide the contents of the email.
                response = client.send_email(
                    Destination={
                        'ToAddresses': [
                            RECIPIENT,
                        ],
                    },
                    Message={
                        'Body': {
                            'Html': {
                                'Charset': CHARSET,
                                'Data': FORMATTED_HTML,
                            },
                            'Text': {
                                'Charset': CHARSET,
                                'Data': BODY_TEXT,
                            },
                        },
                        'Subject': {
                            'Charset': CHARSET,
                            'Data': SUBJECT,
                        }
                    },
                    Source=SENDER,
                )
            # Display an error if something goes wrong.
            except ClientError as e:
                logger.warn("ClientError while notifying.", exc_info=1, client_error_message=e.response['Error']['Message'])
                job_context['job'].success = False
                job_context['job'].failure_reason = e.response['Error']['Message']
                job_context['success'] = False
                return job_context
            except Exception as e:
                logger.warn("General failure when trying to send email.", exc_info=1, result_url=job_context["result_url"])
                job_context['job'].success = False
                job_context['job'].failure_reason = str(e)
                job_context['success'] = False
                return job_context

            job_context["dataset"].email_sent = True
            job_context["dataset"].save()

    return job_context

def _update_result_objects(job_context: Dict) -> Dict:
    """Closes out the dataset object."""

    dataset = job_context["dataset"]
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.is_available = True
    dataset.expires_on = timezone.now() + timedelta(days=7)
    dataset.save()

    job_context['success'] = True

    return job_context

def smash(job_id: int, upload=True) -> None:
    """ Main Smasher interface """

    pipeline = Pipeline(name=PipelineEnum.SMASHER.value)
    return utils.run_pipeline({ "job_id": job_id,
                                "upload": upload,
                                "pipeline": pipeline
                            },
                       [utils.start_job,
                        smashing_utils.prepare_files,
                        _smash_all,
                        _upload,
                        _notify,
                        _update_result_objects,
                        utils.end_job])

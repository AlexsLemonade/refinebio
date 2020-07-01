# -*- coding: utf-8 -*-

import logging
import os
import shutil
import time
from datetime import timedelta
from pathlib import Path
from typing import Dict, List
from urllib.parse import quote

from django.conf import settings
from django.utils import timezone

import boto3
import pandas as pd
import psutil
import requests
from botocore.exceptions import ClientError
from sklearn import preprocessing

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    OriginalFile,
    Pipeline,
    SampleResultAssociation,
)
from data_refinery_common.utils import calculate_file_size, calculate_sha1, get_env_variable
from data_refinery_workers.processors import smashing_utils, utils

RESULTS_BUCKET = get_env_variable("S3_RESULTS_BUCKET_NAME", "refinebio-results-bucket")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
BODY_HTML = (
    Path("data_refinery_workers/processors/smasher_email.min.html").read_text().replace("\n", "")
)
BODY_ERROR_HTML = (
    Path("data_refinery_workers/processors/smasher_email_error.min.html")
    .read_text()
    .replace("\n", "")
)
BYTES_IN_GB = 1024 * 1024 * 1024
logger = get_and_configure_logger(__name__)
### DEBUG ###
logger.setLevel(logging.getLevelName("DEBUG"))

PROCESS_POOL_SIZE = max(1, int(psutil.cpu_count() / 2 - 1))

SCALERS = {
    "MINMAX": preprocessing.MinMaxScaler,
    "STANDARD": preprocessing.StandardScaler,
    "ROBUST": preprocessing.RobustScaler,
}


def log_state(message, job_id, start_time=False):
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process(os.getpid())
        ram_in_GB = process.memory_info().rss / BYTES_IN_GB
        logger.debug(message, total_cpu=psutil.cpu_percent(), process_ram=ram_in_GB, job_id=job_id)

        if start_time:
            logger.debug("Duration: %s" % (time.time() - start_time), job_id=job_id)
        else:
            return time.time()


def _inner_join(job_context: Dict) -> pd.DataFrame:
    """Performs an inner join across the all_frames key of job_context.

    Returns a dataframe, not the job_context.

    TODO: This function should be mostly unnecessary now because we
    pretty much do this in the smashing utils but I don't want to rip
    it out right now .
    """
    # Merge all of the frames we've gathered into a single big frame, skipping duplicates.
    # TODO: If the very first frame is the wrong platform, are we boned?
    merged = job_context["all_frames"][0]
    i = 1

    old_len_merged = len(merged)
    merged_backup = merged

    while i < len(job_context["all_frames"]):
        frame = job_context["all_frames"][i]
        i = i + 1

        if i % 1000 == 0:
            logger.info("Smashing keyframe", i=i, job_id=job_context["job"].id)

        # I'm not sure where these are sneaking in from, but we don't want them.
        # Related: https://github.com/AlexsLemonade/refinebio/issues/390
        breaker = False
        for column in frame.columns:
            if column in merged.columns:
                breaker = True

        if breaker:
            logger.warning(
                "Column repeated for smash job!",
                dataset_id=job_context["dataset"].id,
                job_id=job_context["job"].id,
                column=column,
            )
            continue

        # This is the inner join, the main "Smash" operation
        merged = merged.merge(frame, how="inner", left_index=True, right_index=True)

        new_len_merged = len(merged)
        if new_len_merged < old_len_merged:
            logger.warning(
                "Dropped rows while smashing!",
                dataset_id=job_context["dataset"].id,
                old_len_merged=old_len_merged,
                new_len_merged=new_len_merged,
            )
        if new_len_merged == 0:
            logger.warning(
                "Skipping a bad merge frame!",
                dataset_id=job_context["dataset"].id,
                job_id=job_context["job"].id,
                old_len_merged=old_len_merged,
                new_len_merged=new_len_merged,
                bad_frame_number=i,
            )
            merged = merged_backup
            new_len_merged = len(merged)
            try:
                job_context["unsmashable_files"].append(frame.columns[0])
            except Exception:
                # Something is really, really wrong with this frame.
                pass

        old_len_merged = len(merged)
        merged_backup = merged

    return merged


def process_frames_for_key(key: str, input_files: List[ComputedFile], job_context: Dict) -> Dict:
    """Download, read, and chunk processed sample files from s3.

    `key` is the species or experiment whose samples are contained in `input_files`.

    Will add to job_context the key 'all_frames', a list of pandas
    dataframes containing all the samples' data. Also adds the key
    'unsmashable_files' containing a list of paths that were
    determined to be unsmashable.
    """
    job_context["original_merged"] = pd.DataFrame()

    start_all_frames = log_state(
        "Building list of all_frames key {}".format(key), job_context["job"].id
    )

    job_context["all_frames"] = []
    for (computed_file, sample) in input_files:
        frame_data = smashing_utils.process_frame(
            job_context["work_dir"],
            computed_file,
            sample.accession_code,
            job_context["dataset"].aggregate_by,
        )

        if frame_data is not None:
            job_context["all_frames"].append(frame_data)
        else:
            logger.warning(
                "Unable to smash file",
                computed_file=computed_file.id,
                dataset_id=job_context["dataset"].id,
                job_id=job_context["job"].id,
            )
            job_context["unsmashable_files"].append(computed_file.filename)

    log_state(
        "Finished building list of all_frames key {}".format(key),
        job_context["job"].id,
        start_all_frames,
    )

    return job_context


def _smash_key(job_context: Dict, key: str, input_files: List[ComputedFile]) -> Dict:
    """Smash all of the input files together for a given key.

    Steps:
        Combine common genes (pandas merge)
        Transpose such that genes are columns (features)
        Scale features with sci-kit learn
        Transpose again such that samples are columns and genes are rows
    """
    start_smash = log_state("start _smash_key for {}".format(key), job_context["job"].id)

    # Check if we need to copy the quant.sf files
    if job_context["dataset"].quant_sf_only:
        outfile_dir = job_context["output_dir"] + key + "/"
        os.makedirs(outfile_dir, exist_ok=True)
        samples = [sample for (_, sample) in input_files]
        job_context["num_samples"] += smashing_utils.sync_quant_files(
            outfile_dir, samples, job_context["filtered_samples"]
        )
        # we ONLY want to give quant sf files to the user if that's what they requested
        return job_context

    job_context = process_frames_for_key(key, input_files, job_context)

    if len(job_context["all_frames"]) < 1:
        logger.error(
            "Was told to smash a key with no frames!", job_id=job_context["job"].id, key=key
        )
        # TODO: is this the proper way to handle this? I can see us
        # not wanting to fail an entire dataset because one experiment
        # had a problem, but I also think it could be problematic to
        # just skip an experiment and pretend nothing went wrong.
        return job_context

    merged = _inner_join(job_context)

    job_context["original_merged"] = merged
    log_state("end build all frames", job_context["job"].id, start_smash)
    start_qn = log_state("start qn", job_context["job"].id, start_smash)

    # Quantile Normalization
    if job_context["dataset"].quantile_normalize:
        job_context["merged_no_qn"] = merged
        job_context["organism"] = job_context["dataset"].get_samples()[0].organism
        job_context = smashing_utils.quantile_normalize(job_context)
        merged = job_context.get("merged_qn", None)

    # End QN
    log_state("end qn", job_context["job"].id, start_qn)

    # Transpose before scaling
    transposed = merged.transpose()

    start_scaler = log_state("starting scaler", job_context["job"].id)
    # Scaler
    if job_context["dataset"].scale_by != "NONE":
        scale_funtion = SCALERS[job_context["dataset"].scale_by]
        scaler = scale_funtion(copy=True)
        scaler.fit(transposed)
        scaled = pd.DataFrame(
            scaler.transform(transposed), index=transposed.index, columns=transposed.columns
        )
        # Untranspose
        untransposed = scaled.transpose()
    else:
        # Wheeeeeeeeeee
        untransposed = transposed.transpose()
    log_state("end scaler", job_context["job"].id, start_scaler)

    # This is just for quality assurance in tests.
    job_context["final_frame"] = untransposed

    # Normalize the Header format
    untransposed.index.rename("Gene", inplace=True)

    outfile_dir = job_context["output_dir"] + key + "/"
    os.makedirs(outfile_dir, exist_ok=True)
    outfile = outfile_dir + key + ".tsv"
    job_context["smash_outfile"] = outfile
    untransposed.to_csv(outfile, sep="\t", encoding="utf-8")

    log_state("end _smash_key for {}".format(key), job_context["job"].id, start_smash)

    return job_context


def _smash_all(job_context: Dict) -> Dict:
    """Perform smashing on all species/experiments in the dataset.
    """
    start_smash = log_state("start smash", job_context["job"].id)

    job_context["unsmashable_files"] = []
    job_context["num_samples"] = 0

    # Smash all of the sample sets
    logger.debug(
        "About to smash!",
        dataset_count=len(job_context["dataset"].data),
        job_id=job_context["job"].id,
    )

    try:
        # Once again, `key` is either a species name or an experiment accession
        for key, input_files in job_context.pop("input_files").items():
            job_context = _smash_key(job_context, key, input_files)
    except Exception as e:
        raise utils.ProcessorJobError(
            "Could not smash dataset: " + str(e),
            success=False,
            dataset_id=job_context["dataset"].id,
            num_input_files=job_context["num_input_files"],
        )

    smashing_utils.write_non_data_files(job_context)

    # Finally, compress all files into a zip
    final_zip_base = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk)
    try:
        shutil.make_archive(final_zip_base, "zip", job_context["output_dir"])
    except:
        raise utils.ProcessorJobError("Smash Error while generating zip file", success=False)

    job_context["output_file"] = final_zip_base + ".zip"

    job_context["dataset"].success = True
    job_context["dataset"].save()

    logger.debug("Created smash output!", archive_location=job_context["output_file"])

    log_state("end smash", job_context["job"].id, start_smash)
    return job_context


def _upload(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """
    if not job_context.get("upload", True) or not settings.RUNNING_IN_CLOUD:
        return job_context

    s3_client = boto3.client("s3")
    output_filename = job_context["output_file"].split("/")[-1]

    try:
        # Note that file expiry is handled by the S3 object lifecycle,
        # managed by terraform.
        s3_client.upload_file(
            job_context["output_file"],
            RESULTS_BUCKET,
            output_filename,
            ExtraArgs={"ACL": "public-read"},
        )
    except Exception:
        raise utils.ProcessorJobError(
            "Failed to upload smash result file.", success=False, file=job_context["output_file"]
        )

    result_url = "https://s3.amazonaws.com/" + RESULTS_BUCKET + "/" + output_filename

    job_context["result_url"] = result_url

    logger.debug("Result uploaded!", result_url=job_context["result_url"])

    return job_context


def _notify(job_context: Dict) -> Dict:
    """ Use AWS SES to notify a user of a smash result.. """

    if not job_context.get("upload", True) or not settings.RUNNING_IN_CLOUD:
        return job_context

    # Send a notification to slack when a dataset fails to be processed
    if job_context["job"].success is False:
        try:
            _notify_slack_failed_dataset(job_context)
        except Exception as e:
            logger.warn(e)  # It doesn't really matter if this didn't work

    # Don't send an email if we don't have address.
    if job_context["dataset"].email_address:
        # Try to send the email.
        try:
            _notify_send_email(job_context)
        # Display an error if something goes wrong.
        except ClientError as e:
            raise utils.ProcessorJobError(
                "ClientError while notifying",
                success=False,
                exc_info=1,
                client_error_message=e.response["Error"]["Message"],
            )
        except Exception:
            raise utils.ProcessorJobError(
                "General failure when trying to send email.",
                success=False,
                exc_info=1,
                result_url=job_context["result_url"],
            )

    # We don't want to retry this dataset after we send a notification to users
    # https://github.com/alexslemonade/refinebio/issues/1944
    job_context["job"].no_retry = True
    job_context["job"].save()

    return job_context


def _notify_slack_failed_dataset(job_context: Dict):
    """ Send a slack notification when a dataset fails to smash """

    # Link to the dataset page, where the user can re-try the download job
    dataset_url = "https://www.refine.bio/dataset/" + str(job_context["dataset"].id)

    requests.post(
        settings.ENGAGEMENTBOT_WEBHOOK,
        json={
            "username": "EngagementBot",
            "icon_emoji": ":halal:",
            "attachments": [
                {
                    "fallback": "Dataset failed processing.",
                    "title": "Dataset failed processing",
                    "title_link": dataset_url,
                    "color": "#db3b28",
                    "text": job_context["job"].failure_reason,
                    "fields": [
                        {
                            "title": "Dataset id",
                            "value": str(job_context["dataset"].id),
                            "short": True,
                        },
                        {
                            "title": "Email",
                            "value": job_context["dataset"].email_address,
                            "short": True,
                        },
                    ],
                    "footer": "Refine.bio",
                    "footer_icon": "https://s3.amazonaws.com/refinebio-email/logo-2x.png",
                }
            ],
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )


def _notify_send_email(job_context):
    """ Send email notification to the user if the dataset succeded or failed. """
    dataset_url = "https://www.refine.bio/dataset/" + str(job_context["dataset"].id)

    SENDER = "Refine.bio Mail Robot <noreply@refine.bio>"
    RECIPIENT = job_context["dataset"].email_address
    AWS_REGION = "us-east-1"
    CHARSET = "UTF-8"

    if job_context["job"].success is False:
        SUBJECT = "There was a problem processing your refine.bio dataset :("
        BODY_TEXT = (
            "We tried but were unable to process your requested dataset. Error was: \n\n"
            + str(job_context["job"].failure_reason)
            + "\nDataset ID: "
            + str(job_context["dataset"].id)
            + "\n We have been notified and are looking into the problem. \n\nSorry!"
        )

        ERROR_EMAIL_TITLE = quote("I can't download my dataset")
        ERROR_EMAIL_BODY = quote(
            """
        [What browser are you using?]
        [Add details of the issue you are facing]

        ---
        """
            + str(job_context["dataset"].id)
        )

        FORMATTED_HTML = (
            BODY_ERROR_HTML.replace("REPLACE_DATASET_URL", dataset_url)
            .replace("REPLACE_ERROR_TEXT", job_context["job"].failure_reason)
            .replace(
                "REPLACE_NEW_ISSUE",
                "https://github.com/AlexsLemonade/refinebio/issues/new?title={0}&body={1}&labels=bug".format(
                    ERROR_EMAIL_TITLE, ERROR_EMAIL_BODY
                ),
            )
            .replace(
                "REPLACE_MAILTO",
                "mailto:ccdl@alexslemonade.org?subject={0}&body={1}".format(
                    ERROR_EMAIL_TITLE, ERROR_EMAIL_BODY
                ),
            )
        )
        job_context["success"] = False
    else:
        SUBJECT = "Your refine.bio Dataset is Ready!"
        BODY_TEXT = "Hot off the presses:\n\n" + dataset_url + "\n\nLove!,\nThe refine.bio Team"
        FORMATTED_HTML = BODY_HTML.replace("REPLACE_DOWNLOAD_URL", dataset_url).replace(
            "REPLACE_DATASET_URL", dataset_url
        )

    # Create a new SES resource and specify a region.
    client = boto3.client("ses", region_name=AWS_REGION)

    # Provide the contents of the email.
    response = client.send_email(
        Destination={"ToAddresses": [RECIPIENT,],},
        Message={
            "Body": {
                "Html": {"Charset": CHARSET, "Data": FORMATTED_HTML,},
                "Text": {"Charset": CHARSET, "Data": BODY_TEXT,},
            },
            "Subject": {"Charset": CHARSET, "Data": SUBJECT,},
        },
        Source=SENDER,
    )


def _update_result_objects(job_context: Dict) -> Dict:
    """Closes out the dataset object."""
    dataset = job_context["dataset"]

    dataset.s3_bucket = RESULTS_BUCKET
    dataset.s3_key = job_context["output_file"].split("/")[-1]
    dataset.size_in_bytes = calculate_file_size(job_context["output_file"])
    dataset.sha1 = calculate_sha1(job_context["output_file"])
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.is_available = True
    dataset.expires_on = timezone.now() + timedelta(days=7)
    dataset.save()

    if settings.RUNNING_IN_CLOUD and job_context.get("upload", True):
        # File is uploaded and the metadata is updated, can delete the local.
        try:
            os.remove(job_context["output_file"])
        except OSError:
            pass

    job_context["success"] = True

    return job_context


def smash(job_id: int, upload=True) -> None:
    """ Main Smasher interface """
    pipeline = Pipeline(name=PipelineEnum.SMASHER.value)
    job_context = utils.run_pipeline(
        {"job_id": job_id, "upload": upload, "pipeline": pipeline},
        [
            utils.start_job,
            smashing_utils.prepare_files,
            _smash_all,
            _upload,
            _update_result_objects,
            utils.end_job,
        ],
    )
    # ensure that `notify` is always called so that users get emails in case processing fails or succeeds
    job_context = _notify(job_context)
    return job_context

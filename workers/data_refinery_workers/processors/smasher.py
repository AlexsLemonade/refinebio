from __future__ import absolute_import, unicode_literals

import boto3
import json
import os
import shutil
import string
import warnings

from botocore.exceptions import ClientError
from datetime import datetime, timedelta
from django.utils import timezone
from typing import Dict

import pandas as pd
from sklearn import preprocessing

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import OriginalFile, ComputationalResult, ComputedFile, SampleResultAssociation
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.utils import get_env_variable

RESULTS_BUCKET = get_env_variable("S3_RESULTS_BUCKET_NAME", "refinebio-results-bucket")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """
    Fetches and prepares the files to smash.
    """

    job_context['input_files'] = {}
    for key, samples in job_context["samples"].items():
        all_sample_files = []
        for sample in samples:
            all_sample_files = all_sample_files + list(sample.get_result_files())
        all_sample_files = list(set(all_sample_files))
        job_context['input_files'][key] = all_sample_files

    return job_context

def _smash(job_context: Dict) -> Dict:
    """
    Smash all of the samples together!

    Steps:
        Combine common genes (pandas merge)
        Transpose such that genes are columns (features)
        Scale features with sci-kit learn
        Transpose again such that samples are columns and genes are rows
    """

    try:
        # Prepare the output directory
        smash_path = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk) + "/"
        os.makedirs(smash_path, exist_ok=True)

        scalers = {
            'MINMAX': preprocessing.MinMaxScaler,
            'STANDARD': preprocessing.StandardScaler,
            'ROBUST': preprocessing.RobustScaler,
        }

        num_samples = 0
        # Smash all of the sample sets
        for key, input_files in job_context['input_files'].items():

            print(key)
            print(input_files)

            # Merge all the frames into one
            all_frames = []
            for computed_file in input_files:

                print(computed_file.absolute_file_path)
                print(os.path.exists(computed_file.absolute_file_path))

                data = pd.DataFrame.from_csv(str(computed_file.absolute_file_path), sep='\t', header=0)
                all_frames.append(data)
                num_samples = num_samples + 1
            merged = all_frames[0]
            i = 1
            while i < len(all_frames):
                merged = merged.merge(all_frames[i], left_index=True, right_index=True)
                i = i + 1

            # Transpose before scaling
            transposed = merged.transpose()

            # Scaler
            if job_context['dataset'].scale_by != "NONE":
                scale_funtion = scalers[job_context['dataset'].scale_by]
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

            # Write to temp file with dataset UUID in filename.
            outfile = smash_path + key + ".tsv"
            untransposed.to_csv(outfile, sep='\t', encoding='utf-8')

        # Copy LICENSE.txt and README.md files
        shutil.copy("README_DATASET.md", smash_path + "README.md")
        shutil.copy("LICENSE_DATASET.txt", smash_path + "LICENSE.TXT")

        # Create metadata file.
        metadata = {}
        metadata['files'] = os.listdir(smash_path)
        metadata['num_samples'] = num_samples
        metadata['num_experiments'] = job_context["experiments"].count()
        metadata['aggregate_by'] = job_context["dataset"].aggregate_by
        metadata['scale_by'] = job_context["dataset"].scale_by

        samples = {}
        for sample in job_context["dataset"].get_samples():
            samples[sample.title] = sample.to_metadata_dict()
        metadata['samples'] = samples

        experiments = {}
        for experiment in job_context["dataset"].get_experiments():
            experiments[experiment.accession_code] = experiment.to_metadata_dict()
        metadata['experiments'] = experiments

        metadata['created_at'] = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
        with open(smash_path + 'metadata.json', 'w') as metadata_file: 
            json.dump(metadata, metadata_file, indent=4, sort_keys=True)

        # Finally, compress all files into a zip
        final_zip = smash_path + str(job_context['dataset'].id)
        shutil.make_archive(final_zip, 'zip', smash_path)
        job_context["output_file"] = final_zip + ".zip"
    except Exception as e:
        job_context['dataset'].success = False
        job_context['dataset'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].save()
        job_context['success'] = False
        job_context['failure_reason'] = str(e)
        return job_context

    job_context['dataset'].success = True
    job_context['dataset'].save()

    return job_context

def _upload_and_notify(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """

    if job_context.get("upload", True):

        ##
        # S3
        ##
        s3_client = boto3.client('s3')

        # Note that file expiry is handled by the S3 object lifecycle,
        # managed by terraform.
        s3_client.upload_file(
                job_context["output_file"], 
                RESULTS_BUCKET, 
                job_context["output_file"].split('/')[-1],
                ExtraArgs={'ACL':'public-read'}
            )
        result_url = "https://s3.amazonaws.com/" + RESULTS_BUCKET + "/" + job_context["output_file"].split('/')[-1]
        job_context["result_url"] = result_url

        job_context["dataset"].s3_bucket = RESULTS_BUCKET
        job_context["dataset"].s3_key = job_context["output_file"].split('/')[-1]
        job_context["dataset"].save()

        ##
        # SES
        ##

        # Don't send an email if we don't have address.
        if job_context["dataset"].email_address:
            # XXX: This address must be verified with Amazon SES.
            SENDER = "Refine.bio Mail Robot <noreply@refine.bio>"
            RECIPIENT = job_context["dataset"].email_address
            AWS_REGION = "us-east-1"
            SUBJECT = "Your refine.bio Dataset is Ready!"
            BODY_TEXT = "Hot off the presses:\n\n" + result_url + "\n\nLove!,\nThe refine.bio Team"
            BODY_HTML = "Hot off the presses:<br /><br />" + result_url + "<br /><br />Love!,<br />The refine.bio Team"
            CHARSET = "UTF-8"

            # Create a new SES resource and specify a region.
            client = boto3.client('ses', region_name=AWS_REGION)

            # Try to send the email.
            try:
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
                                'Data': BODY_HTML,
                            },
                            'Text': {
                                'Charset': CHARSET,
                                'Data': BODY_TEXT,
                            },
                        },
                        'Subject': {
                            'Charset': CHARSET,
                            'Data': SUBJECT,
                        },
                    },
                    Source=SENDER,
                )
            # Display an error if something goes wrong. 
            except ClientError as e:
                logger.error(e.response['Error']['Message'])
            else:
                job_context["dataset"].email_sent = True
                job_context["dataset"].save()

    return job_context

def _update_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a run is complete """

    dataset = job_context["dataset"]
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.is_available = True
    dataset.expires_on = timezone.now() + timedelta(days=1)
    dataset.save()

    return job_context

def smash(job_id: int, upload=True) -> None:
    """ Main Smasher interface """

    return utils.run_pipeline({"job_id": job_id, "upload": upload},
                       [utils.start_job,
                        _prepare_files,
                        _smash,
                        _upload_and_notify,
                        _update_result_objects,
                        utils.end_job])

from __future__ import absolute_import, unicode_literals

import boto3
import csv
import simplejson as json
import os
import shutil
import string
import warnings

from botocore.exceptions import ClientError
from datetime import datetime, timedelta
from django.utils import timezone
from django.conf import settings
from typing import Dict

import numpy as np
import pandas as pd
from sklearn import preprocessing

import rpy2
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    OriginalFile,
    Pipeline,
    ComputationalResult,
    ComputedFile,
    SampleResultAssociation
)
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

    all_sample_files = []
    job_context['input_files'] = {}

    for key, samples in job_context["samples"].items():
        samples_for_key = []
        for sample in samples:
            samples_for_key = samples_for_key + list(sample.get_result_files())
        samples_for_key = list(set(samples_for_key))
        job_context['input_files'][key] = samples_for_key
        all_sample_files = all_sample_files + samples_for_key

    if all_sample_files == []:
        logger.error("Couldn't get any files to smash for Smash job!!",
            dataset_id=job_context['dataset'].id,
            samples=job_context["samples"])
        job_context['success'] = False
        return job_context

    # So these get deleted from disk after..
    for computed_file in all_sample_files:
        job_context['computed_files'].append(computed_file)

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
        # Ensure we have a fresh smash directory
        shutil.rmtree(smash_path, ignore_errors=True)
        os.makedirs(smash_path, exist_ok=True)

        scalers = {
            'MINMAX': preprocessing.MinMaxScaler,
            'STANDARD': preprocessing.StandardScaler,
            'ROBUST': preprocessing.RobustScaler,
        }

        unsmashable_files = []
        num_samples = 0

        # Smash all of the sample sets
        logger.info("About to smash!",
                input_files=job_context['input_files'],
                dataset_data=job_context['dataset'].data,
            )
        for key, input_files in job_context['input_files'].items():

            # Merge all the frames into one
            all_frames = []

            for computed_file in input_files:

                computed_file_path = str(computed_file.get_synced_file_path())

                # Bail appropriately if this isn't a real file.
                if not os.path.exists(computed_file.get_synced_file_path()):
                    raise ValueError("Smasher received non-existent file path.")

                try:
                    data = pd.read_csv(computed_file_path, sep='\t', header=0, index_col=0, error_bad_lines=False)

                    # Strip any funky whitespace
                    data.columns = data.columns.str.strip()
                    data = data.dropna(axis='columns', how='all')

                    # Make sure the index type is correct
                    data.index = data.index.map(str)

                    # via https://github.com/AlexsLemonade/refinebio/issues/330:
                    #   aggregating by experiment -> return untransformed output from tximport
                    #   aggregating by species -> log2(x + 1) tximport output
                    if job_context['dataset'].aggregate_by == 'SPECIES':
                        if 'lengthScaledTPM' in computed_file_path:
                            data = data + 1
                            data = np.log2(data)

                    # Detect if this data hasn't been log2 scaled yet.
                    # Ideally done in the NO-OPPER, but sanity check here.
                    if ("lengthScaledTPM" not in computed_file_path) and (data.max() > 100).any():
                        logger.info("Detected non-log2 microarray data.", file=computed_file)
                        data = np.log2(data)

                    # Ensure that we don't have any dangling Brainarray-generated probe symbols.
                    # BA likes to leave '_at', signifying probe identifiers,
                    # on their converted, non-probe identifiers. It makes no sense.
                    # So, we chop them off and don't worry about it.
                    data.index = data.index.str.replace('_at', '')

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

                    # Explicitly title this dataframe
                    try:
                        data.columns = [computed_file.samples.all()[0].title]
                    except ValueError:
                        # This sample have have multiple channels, or something else.
                        # Don't mess with it.
                        pass
                    except Exception as e:
                        # Okay, somebody probably forgot to create a SampleComputedFileAssociation
                        data.columns = [computed_file.filename]

                    all_frames.append(data)
                    num_samples = num_samples + 1

                except Exception as e:
                    unsmashable_files.append(computed_file_path)
                    logger.exception("Unable to smash file",
                        file=computed_file_path,
                        dataset_id=job_context['dataset'].id,
                        )
                    continue

            job_context['all_frames'] = all_frames

            if len(all_frames) < 1:
                logger.warning("Was told to smash a frame with no frames!",
                    key=key,
                    input_files=str(input_files)
                )
                continue

            # Merge all of the frames we've gathered into a single big frame, skipping duplicates.
            merged = all_frames[0]
            i = 1
            while i < len(all_frames):
                frame = all_frames[i]
                i = i + 1

                # I'm not sure where these are sneaking in from, but we don't want them.
                # Related: https://github.com/AlexsLemonade/refinebio/issues/390
                breaker = False
                for column in frame.columns:
                    if column in merged.columns:
                        breaker = True
                if breaker:
                    continue

                # This is the inner join, the main "Smash" operation
                merged = merged.merge(frame, left_index=True, right_index=True)

            job_context['original_merged'] = merged

            # Quantile Normalization
            if job_context['dataset'].quantile_normalize:
                try:
                    # Prepare our QN target file
                    organism = computed_file.samples.first().organism
                    qn_target = utils.get_most_recent_qn_target_for_organism(organism)

                    if not qn_target:
                        logger.error("Could not find QN target for Organism!",
                            organism=organism,
                            dataset_id=job_context['dataset'].id,
                            dataset_data=job_context['dataset'].data,
                        )
                    else:
                        qn_target_path = qn_target.sync_from_s3()
                        qn_target_frame = pd.read_csv(qn_target_path, sep='\t', header=None, index_col=None, error_bad_lines=False)

                        # Prepare our RPy2 bridge
                        pandas2ri.activate()
                        preprocessCore = importr('preprocessCore')
                        as_numeric = rlang("as.numeric")
                        data_matrix = rlang('data.matrix')

                        # Convert the smashed frames to an R numeric Matrix
                        # and the target Dataframe into an R numeric Vector
                        target_vector = as_numeric(qn_target_frame[0])
                        merged_matrix = data_matrix(merged)

                        # Perform the Actual QN
                        reso = preprocessCore.normalize_quantiles_use_target(
                                                            x=merged_matrix,
                                                            target=target_vector,
                                                            copy=True
                                                        )

                        # And finally convert back to Pandas
                        ar = np.array(reso)
                        new_merged = pd.DataFrame(ar, columns=merged.columns, index=merged.index)
                        job_context['merged_no_qn'] = merged
                        job_context['merged_qn'] = new_merged
                        merged = new_merged
                except Exception as e:
                    logger.exception("Problem occured during quantile normalization",
                        dataset_id=job_context['dataset'].id,
                        dataset_data=job_context['dataset'].data,
                    )

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

            # This is just for quality assurance in tests.
            job_context['final_frame'] = untransposed

            # Write to temp file with dataset UUID in filename.
            outfile = smash_path + key + ".tsv"
            untransposed.to_csv(outfile, sep='\t', encoding='utf-8')

        # Copy LICENSE.txt and README.md files
        shutil.copy("README_DATASET.md", smash_path + "README.md")
        shutil.copy("LICENSE_DATASET.txt", smash_path + "LICENSE.TXT")

        # Create metadata file.
        metadata = {}

        metadata['num_samples'] = num_samples
        metadata['num_experiments'] = job_context["experiments"].count()
        metadata['aggregate_by'] = job_context["dataset"].aggregate_by
        metadata['scale_by'] = job_context["dataset"].scale_by
        # https://github.com/AlexsLemonade/refinebio/pull/421#discussion_r203799646
        metadata['non_aggregated_files'] = unsmashable_files

        samples = {}
        for sample in job_context["dataset"].get_samples():
            samples[sample.title] = sample.to_metadata_dict()
        metadata['samples'] = samples

        experiments = {}
        for experiment in job_context["dataset"].get_experiments():
            exp_dict = experiment.to_metadata_dict()
            exp_dict['sample_titles'] = [v for v in experiment.samples.all().values_list('title', flat=True)]
            experiments[experiment.accession_code] = exp_dict
        metadata['experiments'] = experiments

        # Metadata to TSV
        if job_context["dataset"].aggregate_by == "EXPERIMENT":
            for title, keys in metadata['experiments'].items():
                with open(smash_path + title + '_metadata.tsv', 'w') as output_file:
                    sample_titles = keys['sample_titles']
                    keys = list(metadata['samples'][sample_titles[0]].keys()) + ['sample_id']
                    dw = csv.DictWriter(output_file, keys, delimiter='\t')
                    dw.writeheader()
                    for key, value in metadata['samples'].items():
                        if key not in sample_titles:
                            continue
                        else:
                            value['sample_id'] = key
                            dw.writerow(value)
        else:
            with open(smash_path + 'metadata.tsv', 'w') as output_file:
                sample_ids = list(metadata['samples'].keys())
                keys = list(metadata['samples'][sample_ids[0]].keys()) + ['sample_id']
                dw = csv.DictWriter(output_file, keys, delimiter='\t')
                dw.writeheader()
                for key, value in metadata['samples'].items():
                    value['sample_id'] = key
                    dw.writerow(value)

        metadata['files'] = os.listdir(smash_path)

        # Metadata to JSON
        metadata['created_at'] = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
        with open(smash_path + 'metadata.json', 'w') as metadata_file:
            json.dump(metadata, metadata_file, indent=4, sort_keys=True)

        # Finally, compress all files into a zip
        final_zip_base = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk)
        shutil.make_archive(final_zip_base, 'zip', smash_path)
        job_context["output_file"] = final_zip_base + ".zip"
        # and clean up the unzipped directory.
        shutil.rmtree(smash_path)
    except Exception as e:
        logger.exception("Could not smash dataset.",
                        dataset_id=job_context['dataset'].id,
                        job_id=job_context['job_id'],
                        input_files=job_context['input_files'])
        job_context['dataset'].success = False
        job_context['job'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].save()
        # Delay failing this pipeline until the failure notify has been sent
        # job_context['success'] = False
        job_context['failure_reason'] = str(e)
        return job_context

    job_context['metadata'] = metadata
    job_context['dataset'].success = True
    job_context['dataset'].save()

    logger.info("Created smash output!",
        archive_location=job_context["output_file"])

    return job_context

def _upload(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """

    try:
        if job_context.get("upload", True):
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

            logger.info("Result uploaded!",
                    result_url=job_context["result_url"]
                )

            job_context["dataset"].s3_bucket = RESULTS_BUCKET
            job_context["dataset"].s3_key = job_context["output_file"].split('/')[-1]
            job_context["dataset"].save()

            # File is uploaded, we can delete the local.
            try:
                os.remove(job_context["output_file"])
            except OSError:
                pass

    except Exception as e:
        logger.exception("Failed to upload smash result file.", file=job_context["output_file"])
        job_context['job'].success = False
        job_context['job'].failure_reason = str(e)
        # Delay failing this pipeline until the failure notify has been sent
        # job_context['success'] = False

    return job_context

def _notify(job_context: Dict) -> Dict:
    """ Use AWS SES to notify a user of a smash result.. """

    ##
    # SES
    ##
    if job_context.get("upload", True) and settings.RUNNING_IN_CLOUD:

        # Don't send an email if we don't have address.
        if job_context["dataset"].email_address:
            # XXX: This address must be verified with Amazon SES.
            SENDER = "Refine.bio Mail Robot <noreply@refine.bio>"
            RECIPIENT = job_context["dataset"].email_address
            AWS_REGION = "us-east-1"
            CHARSET = "UTF-8"

            if job_context['job'].failure_reason not in ['', None]:
                SUBJECT = "Your refine.bio Dataset is Ready!"
                BODY_TEXT = "Hot off the presses:\n\n" + job_context["result_url"] + "\n\nLove!,\nThe refine.bio Team"
                BODY_HTML = "Hot off the presses:<br /><br />" + job_context["result_url"] + "<br /><br />Love!,<br />The refine.bio Team"
            else:
                SUBJECT = "There was a problem processing your refine.bio dataset :("
                BODY_TEXT = "We tried but were unable to process your requested dataset. Error was: \n\n" + str(job_context['job'].failure_reason) + "\nDataset ID: " + str(dataset.id) + "\n We have been notified and are looking into the problem. \n\nSorry!"
                BODY_HTML = BODY_TEXT = "We tried but were unable to process your requested dataset. Error was: <br /><br />" + job_context['job'].failure_reason + "<br />Dataset: " + str(dataset.id) + "<br /> We have been notified and are looking into the problem. <br /><br />Sorry!"
                job_context['success'] = False

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
                        }
                    },
                    Source=SENDER,
                )
            # Display an error if something goes wrong.
            except ClientError as e:
                logger.exception("ClientError while notifying.", client_error_message=e.response['Error']['Message'])
                job_context['job'].success = False
                job_context['job'].failure_reason = e.response['Error']['Message']
                job_context['success'] = False
                return job_context
            except Exception as e:
                logger.exception("General failure when trying to send email.", result_url=job_context["result_url"])
                job_context['job'].success = False
                job_context['job'].failure_reason = str(e)
                job_context['success'] = False
                return job_context

            job_context["dataset"].email_sent = True
            job_context["dataset"].save()

    # Handle non-cloud too
    if job_context['job'].failure_reason:
        job_context['success'] = False

    return job_context

def _update_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a run is complete """

    dataset = job_context["dataset"]
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.is_available = True
    dataset.expires_on = timezone.now() + timedelta(days=1)
    dataset.save()

    job_context['success'] = True

    return job_context

def _delete_local_files(job_context: Dict) -> Dict:
    """ Removes all of the ComputedFiles that have been synced from S3 """
    for key, input_files in job_context['input_files'].items():
        for computed_file in input_files:
            computed_file.delete_local_file()

    job_context['success'] = True
    return job_context

def smash(job_id: int, upload=True) -> None:
    """ Main Smasher interface """

    pipeline = Pipeline(name=utils.PipelineEnum.SMASHER.value)
    return utils.run_pipeline({ "job_id": job_id, 
                                "upload": upload, 
                                "pipeline": pipeline
                            },
                       [utils.start_job,
                        _prepare_files,
                        _smash,
                        _upload,
                        _notify,
                        _update_result_objects,
                        _delete_local_files,
                        utils.end_job])


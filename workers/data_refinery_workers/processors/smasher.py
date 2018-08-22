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

    job_context['input_files'] = {}

    for key, samples in job_context["samples"].items():
        all_sample_files = []
        for sample in samples:
            all_sample_files = all_sample_files + list(sample.get_result_files())
        all_sample_files = list(set(all_sample_files))
        job_context['input_files'][key] = all_sample_files

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

                merged = merged.merge(frame, left_index=True, right_index=True)

            job_context['merged'] = merged

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
        job_context['success'] = False
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
        job_context['success'] = False

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
            SUBJECT = "Your refine.bio Dataset is Ready!"
            BODY_TEXT = "Hot off the presses:\n\n" + job_context["result_url"] + "\n\nLove!,\nThe refine.bio Team"
            FORMATTED_HTML = BODY_HTML.replace('REPLACEME', job_context["result_url"])
            CHARSET = "UTF-8"

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
    return utils.run_pipeline({"job_id": job_id, "upload": upload, "pipeline": pipeline},
                       [utils.start_job,
                        _prepare_files,
                        _smash,
                        _upload,
                        _notify,
                        _update_result_objects,
                        _delete_local_files,
                        utils.end_job])

# Look ma, I'm React!
BODY_HTML = """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <meta name="viewport" content="width=device-width">
    <title>Title</title>
    <style>
      .wrapper {
  width: 100%; }

#outlook a {
  padding: 0; }

body {
  width: 100% !important;
  min-width: 100%;
  -webkit-text-size-adjust: 100%;
  -ms-text-size-adjust: 100%;
  margin: 0;
  Margin: 0;
  padding: 0;
  -moz-box-sizing: border-box;
  -webkit-box-sizing: border-box;
  box-sizing: border-box; }

.ExternalClass {
  width: 100%; }
  .ExternalClass,
  .ExternalClass p,
  .ExternalClass span,
  .ExternalClass font,
  .ExternalClass td,
  .ExternalClass div {
    line-height: 100%; }

#backgroundTable {
  margin: 0;
  Margin: 0;
  padding: 0;
  width: 100% !important;
  line-height: 100% !important; }

img {
  outline: none;
  text-decoration: none;
  -ms-interpolation-mode: bicubic;
  width: auto;
  max-width: 100%;
  clear: both;
  display: block; }

center {
  width: 100%;
  min-width: 580px; }

a img {
  border: none; }

p {
  margin: 0 0 0 10px;
  Margin: 0 0 0 10px; }

table {
  border-spacing: 0;
  border-collapse: collapse; }

td {
  word-wrap: break-word;
  -webkit-hyphens: auto;
  -moz-hyphens: auto;
  hyphens: auto;
  border-collapse: collapse !important; }

table, tr, td {
  padding: 0;
  vertical-align: top;
  text-align: left; }

@media only screen {
  html {
    min-height: 100%;
    background: #f3f3f3; } }

table.body {
  background: #f3f3f3;
  height: 100%;
  width: 100%; }

table.container {
  background: #fefefe;
  width: 580px;
  margin: 0 auto;
  Margin: 0 auto;
  text-align: inherit; }

table.row {
  padding: 0;
  width: 100%;
  position: relative; }

table.spacer {
  width: 100%; }
  table.spacer td {
    mso-line-height-rule: exactly; }

table.container table.row {
  display: table; }

td.columns,
td.column,
th.columns,
th.column {
  margin: 0 auto;
  Margin: 0 auto;
  padding-left: 16px;
  padding-bottom: 16px; }
  td.columns .column,
  td.columns .columns,
  td.column .column,
  td.column .columns,
  th.columns .column,
  th.columns .columns,
  th.column .column,
  th.column .columns {
    padding-left: 0 !important;
    padding-right: 0 !important; }
    td.columns .column center,
    td.columns .columns center,
    td.column .column center,
    td.column .columns center,
    th.columns .column center,
    th.columns .columns center,
    th.column .column center,
    th.column .columns center {
      min-width: none !important; }

td.columns.last,
td.column.last,
th.columns.last,
th.column.last {
  padding-right: 16px; }

td.columns table:not(.button),
td.column table:not(.button),
th.columns table:not(.button),
th.column table:not(.button) {
  width: 100%; }

td.large-1,
th.large-1 {
  width: 32.33333px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-1.first,
th.large-1.first {
  padding-left: 16px; }

td.large-1.last,
th.large-1.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-1,
.collapse > tbody > tr > th.large-1 {
  padding-right: 0;
  padding-left: 0;
  width: 48.33333px; }

.collapse td.large-1.first,
.collapse th.large-1.first,
.collapse td.large-1.last,
.collapse th.large-1.last {
  width: 56.33333px; }

td.large-1 center,
th.large-1 center {
  min-width: 0.33333px; }

.body .columns td.large-1,
.body .column td.large-1,
.body .columns th.large-1,
.body .column th.large-1 {
  width: 8.33333%; }

td.large-2,
th.large-2 {
  width: 80.66667px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-2.first,
th.large-2.first {
  padding-left: 16px; }

td.large-2.last,
th.large-2.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-2,
.collapse > tbody > tr > th.large-2 {
  padding-right: 0;
  padding-left: 0;
  width: 96.66667px; }

.collapse td.large-2.first,
.collapse th.large-2.first,
.collapse td.large-2.last,
.collapse th.large-2.last {
  width: 104.66667px; }

td.large-2 center,
th.large-2 center {
  min-width: 48.66667px; }

.body .columns td.large-2,
.body .column td.large-2,
.body .columns th.large-2,
.body .column th.large-2 {
  width: 16.66667%; }

td.large-3,
th.large-3 {
  width: 129px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-3.first,
th.large-3.first {
  padding-left: 16px; }

td.large-3.last,
th.large-3.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-3,
.collapse > tbody > tr > th.large-3 {
  padding-right: 0;
  padding-left: 0;
  width: 145px; }

.collapse td.large-3.first,
.collapse th.large-3.first,
.collapse td.large-3.last,
.collapse th.large-3.last {
  width: 153px; }

td.large-3 center,
th.large-3 center {
  min-width: 97px; }

.body .columns td.large-3,
.body .column td.large-3,
.body .columns th.large-3,
.body .column th.large-3 {
  width: 25%; }

td.large-4,
th.large-4 {
  width: 177.33333px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-4.first,
th.large-4.first {
  padding-left: 16px; }

td.large-4.last,
th.large-4.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-4,
.collapse > tbody > tr > th.large-4 {
  padding-right: 0;
  padding-left: 0;
  width: 193.33333px; }

.collapse td.large-4.first,
.collapse th.large-4.first,
.collapse td.large-4.last,
.collapse th.large-4.last {
  width: 201.33333px; }

td.large-4 center,
th.large-4 center {
  min-width: 145.33333px; }

.body .columns td.large-4,
.body .column td.large-4,
.body .columns th.large-4,
.body .column th.large-4 {
  width: 33.33333%; }

td.large-5,
th.large-5 {
  width: 225.66667px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-5.first,
th.large-5.first {
  padding-left: 16px; }

td.large-5.last,
th.large-5.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-5,
.collapse > tbody > tr > th.large-5 {
  padding-right: 0;
  padding-left: 0;
  width: 241.66667px; }

.collapse td.large-5.first,
.collapse th.large-5.first,
.collapse td.large-5.last,
.collapse th.large-5.last {
  width: 249.66667px; }

td.large-5 center,
th.large-5 center {
  min-width: 193.66667px; }

.body .columns td.large-5,
.body .column td.large-5,
.body .columns th.large-5,
.body .column th.large-5 {
  width: 41.66667%; }

td.large-6,
th.large-6 {
  width: 274px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-6.first,
th.large-6.first {
  padding-left: 16px; }

td.large-6.last,
th.large-6.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-6,
.collapse > tbody > tr > th.large-6 {
  padding-right: 0;
  padding-left: 0;
  width: 290px; }

.collapse td.large-6.first,
.collapse th.large-6.first,
.collapse td.large-6.last,
.collapse th.large-6.last {
  width: 298px; }

td.large-6 center,
th.large-6 center {
  min-width: 242px; }

.body .columns td.large-6,
.body .column td.large-6,
.body .columns th.large-6,
.body .column th.large-6 {
  width: 50%; }

td.large-7,
th.large-7 {
  width: 322.33333px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-7.first,
th.large-7.first {
  padding-left: 16px; }

td.large-7.last,
th.large-7.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-7,
.collapse > tbody > tr > th.large-7 {
  padding-right: 0;
  padding-left: 0;
  width: 338.33333px; }

.collapse td.large-7.first,
.collapse th.large-7.first,
.collapse td.large-7.last,
.collapse th.large-7.last {
  width: 346.33333px; }

td.large-7 center,
th.large-7 center {
  min-width: 290.33333px; }

.body .columns td.large-7,
.body .column td.large-7,
.body .columns th.large-7,
.body .column th.large-7 {
  width: 58.33333%; }

td.large-8,
th.large-8 {
  width: 370.66667px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-8.first,
th.large-8.first {
  padding-left: 16px; }

td.large-8.last,
th.large-8.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-8,
.collapse > tbody > tr > th.large-8 {
  padding-right: 0;
  padding-left: 0;
  width: 386.66667px; }

.collapse td.large-8.first,
.collapse th.large-8.first,
.collapse td.large-8.last,
.collapse th.large-8.last {
  width: 394.66667px; }

td.large-8 center,
th.large-8 center {
  min-width: 338.66667px; }

.body .columns td.large-8,
.body .column td.large-8,
.body .columns th.large-8,
.body .column th.large-8 {
  width: 66.66667%; }

td.large-9,
th.large-9 {
  width: 419px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-9.first,
th.large-9.first {
  padding-left: 16px; }

td.large-9.last,
th.large-9.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-9,
.collapse > tbody > tr > th.large-9 {
  padding-right: 0;
  padding-left: 0;
  width: 435px; }

.collapse td.large-9.first,
.collapse th.large-9.first,
.collapse td.large-9.last,
.collapse th.large-9.last {
  width: 443px; }

td.large-9 center,
th.large-9 center {
  min-width: 387px; }

.body .columns td.large-9,
.body .column td.large-9,
.body .columns th.large-9,
.body .column th.large-9 {
  width: 75%; }

td.large-10,
th.large-10 {
  width: 467.33333px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-10.first,
th.large-10.first {
  padding-left: 16px; }

td.large-10.last,
th.large-10.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-10,
.collapse > tbody > tr > th.large-10 {
  padding-right: 0;
  padding-left: 0;
  width: 483.33333px; }

.collapse td.large-10.first,
.collapse th.large-10.first,
.collapse td.large-10.last,
.collapse th.large-10.last {
  width: 491.33333px; }

td.large-10 center,
th.large-10 center {
  min-width: 435.33333px; }

.body .columns td.large-10,
.body .column td.large-10,
.body .columns th.large-10,
.body .column th.large-10 {
  width: 83.33333%; }

td.large-11,
th.large-11 {
  width: 515.66667px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-11.first,
th.large-11.first {
  padding-left: 16px; }

td.large-11.last,
th.large-11.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-11,
.collapse > tbody > tr > th.large-11 {
  padding-right: 0;
  padding-left: 0;
  width: 531.66667px; }

.collapse td.large-11.first,
.collapse th.large-11.first,
.collapse td.large-11.last,
.collapse th.large-11.last {
  width: 539.66667px; }

td.large-11 center,
th.large-11 center {
  min-width: 483.66667px; }

.body .columns td.large-11,
.body .column td.large-11,
.body .columns th.large-11,
.body .column th.large-11 {
  width: 91.66667%; }

td.large-12,
th.large-12 {
  width: 564px;
  padding-left: 8px;
  padding-right: 8px; }

td.large-12.first,
th.large-12.first {
  padding-left: 16px; }

td.large-12.last,
th.large-12.last {
  padding-right: 16px; }

.collapse > tbody > tr > td.large-12,
.collapse > tbody > tr > th.large-12 {
  padding-right: 0;
  padding-left: 0;
  width: 580px; }

.collapse td.large-12.first,
.collapse th.large-12.first,
.collapse td.large-12.last,
.collapse th.large-12.last {
  width: 588px; }

td.large-12 center,
th.large-12 center {
  min-width: 532px; }

.body .columns td.large-12,
.body .column td.large-12,
.body .columns th.large-12,
.body .column th.large-12 {
  width: 100%; }

td.large-offset-1,
td.large-offset-1.first,
td.large-offset-1.last,
th.large-offset-1,
th.large-offset-1.first,
th.large-offset-1.last {
  padding-left: 64.33333px; }

td.large-offset-2,
td.large-offset-2.first,
td.large-offset-2.last,
th.large-offset-2,
th.large-offset-2.first,
th.large-offset-2.last {
  padding-left: 112.66667px; }

td.large-offset-3,
td.large-offset-3.first,
td.large-offset-3.last,
th.large-offset-3,
th.large-offset-3.first,
th.large-offset-3.last {
  padding-left: 161px; }

td.large-offset-4,
td.large-offset-4.first,
td.large-offset-4.last,
th.large-offset-4,
th.large-offset-4.first,
th.large-offset-4.last {
  padding-left: 209.33333px; }

td.large-offset-5,
td.large-offset-5.first,
td.large-offset-5.last,
th.large-offset-5,
th.large-offset-5.first,
th.large-offset-5.last {
  padding-left: 257.66667px; }

td.large-offset-6,
td.large-offset-6.first,
td.large-offset-6.last,
th.large-offset-6,
th.large-offset-6.first,
th.large-offset-6.last {
  padding-left: 306px; }

td.large-offset-7,
td.large-offset-7.first,
td.large-offset-7.last,
th.large-offset-7,
th.large-offset-7.first,
th.large-offset-7.last {
  padding-left: 354.33333px; }

td.large-offset-8,
td.large-offset-8.first,
td.large-offset-8.last,
th.large-offset-8,
th.large-offset-8.first,
th.large-offset-8.last {
  padding-left: 402.66667px; }

td.large-offset-9,
td.large-offset-9.first,
td.large-offset-9.last,
th.large-offset-9,
th.large-offset-9.first,
th.large-offset-9.last {
  padding-left: 451px; }

td.large-offset-10,
td.large-offset-10.first,
td.large-offset-10.last,
th.large-offset-10,
th.large-offset-10.first,
th.large-offset-10.last {
  padding-left: 499.33333px; }

td.large-offset-11,
td.large-offset-11.first,
td.large-offset-11.last,
th.large-offset-11,
th.large-offset-11.first,
th.large-offset-11.last {
  padding-left: 547.66667px; }

td.expander,
th.expander {
  visibility: hidden;
  width: 0;
  padding: 0 !important; }

table.container.radius {
  border-radius: 0;
  border-collapse: separate; }

.block-grid {
  width: 100%;
  max-width: 580px; }
  .block-grid td {
    display: inline-block;
    padding: 8px; }

.up-2 td {
  width: 274px !important; }

.up-3 td {
  width: 177px !important; }

.up-4 td {
  width: 129px !important; }

.up-5 td {
  width: 100px !important; }

.up-6 td {
  width: 80px !important; }

.up-7 td {
  width: 66px !important; }

.up-8 td {
  width: 56px !important; }

table.text-center,
th.text-center,
td.text-center,
h1.text-center,
h2.text-center,
h3.text-center,
h4.text-center,
h5.text-center,
h6.text-center,
p.text-center,
span.text-center {
  text-align: center; }

table.text-left,
th.text-left,
td.text-left,
h1.text-left,
h2.text-left,
h3.text-left,
h4.text-left,
h5.text-left,
h6.text-left,
p.text-left,
span.text-left {
  text-align: left; }

table.text-right,
th.text-right,
td.text-right,
h1.text-right,
h2.text-right,
h3.text-right,
h4.text-right,
h5.text-right,
h6.text-right,
p.text-right,
span.text-right {
  text-align: right; }

span.text-center {
  display: block;
  width: 100%;
  text-align: center; }

@media only screen and (max-width: 596px) {
  .small-float-center {
    margin: 0 auto !important;
    float: none !important;
    text-align: center !important; }
  .small-text-center {
    text-align: center !important; }
  .small-text-left {
    text-align: left !important; }
  .small-text-right {
    text-align: right !important; } }

img.float-left {
  float: left;
  text-align: left; }

img.float-right {
  float: right;
  text-align: right; }

img.float-center,
img.text-center {
  margin: 0 auto;
  Margin: 0 auto;
  float: none;
  text-align: center; }

table.float-center,
td.float-center,
th.float-center {
  margin: 0 auto;
  Margin: 0 auto;
  float: none;
  text-align: center; }

.hide-for-large {
  display: none !important;
  mso-hide: all;
  overflow: hidden;
  max-height: 0;
  font-size: 0;
  width: 0;
  line-height: 0; }
  @media only screen and (max-width: 596px) {
    .hide-for-large {
      display: block !important;
      width: auto !important;
      overflow: visible !important;
      max-height: none !important;
      font-size: inherit !important;
      line-height: inherit !important; } }

table.body table.container .hide-for-large * {
  mso-hide: all; }

@media only screen and (max-width: 596px) {
  table.body table.container .hide-for-large,
  table.body table.container .row.hide-for-large {
    display: table !important;
    width: 100% !important; } }

@media only screen and (max-width: 596px) {
  table.body table.container .callout-inner.hide-for-large {
    display: table-cell !important;
    width: 100% !important; } }

@media only screen and (max-width: 596px) {
  table.body table.container .show-for-large {
    display: none !important;
    width: 0;
    mso-hide: all;
    overflow: hidden; } }

body,
table.body,
h1,
h2,
h3,
h4,
h5,
h6,
p,
td,
th,
a {
  color: #0a0a0a;
  font-family: Helvetica, Arial, sans-serif;
  font-weight: normal;
  padding: 0;
  margin: 0;
  Margin: 0;
  text-align: left;
  line-height: 1.3; }

h1,
h2,
h3,
h4,
h5,
h6 {
  color: inherit;
  word-wrap: normal;
  font-family: Helvetica, Arial, sans-serif;
  font-weight: normal;
  margin-bottom: 10px;
  Margin-bottom: 10px; }

h1 {
  font-size: 34px; }

h2 {
  font-size: 30px; }

h3 {
  font-size: 28px; }

h4 {
  font-size: 24px; }

h5 {
  font-size: 20px; }

h6 {
  font-size: 18px; }

body,
table.body,
p,
td,
th {
  font-size: 16px;
  line-height: 1.3; }

p {
  margin-bottom: 10px;
  Margin-bottom: 10px; }
  p.lead {
    font-size: 20px;
    line-height: 1.6; }
  p.subheader {
    margin-top: 4px;
    margin-bottom: 8px;
    Margin-top: 4px;
    Margin-bottom: 8px;
    font-weight: normal;
    line-height: 1.4;
    color: #8a8a8a; }

small {
  font-size: 80%;
  color: #cacaca; }

a {
  color: #2199e8;
  text-decoration: none;
}
  a:hover {
    color: #147dc2;
  }
  a:active {
    color: #147dc2; }
  a:visited {
    color: #2199e8; }

h1 a,
h1 a:visited,
h2 a,
h2 a:visited,
h3 a,
h3 a:visited,
h4 a,
h4 a:visited,
h5 a,
h5 a:visited,
h6 a,
h6 a:visited {
  color: #2199e8; }

pre {
  background: #f3f3f3;
  margin: 30px 0;
  Margin: 30px 0; }
  pre code {
    color: #cacaca; }
    pre code span.callout {
      color: #8a8a8a;
      font-weight: bold; }
    pre code span.callout-strong {
      color: #ff6908;
      font-weight: bold; }

table.hr {
  width: 100%; }
  table.hr th {
    height: 0;
    max-width: 580px;
    border-top: 0;
    border-right: 0;
    border-bottom: 1px solid #0a0a0a;
    border-left: 0;
    margin: 20px auto;
    Margin: 20px auto;
    clear: both; }

.stat {
  font-size: 40px;
  line-height: 1; }
  p + .stat {
    margin-top: -16px;
    Margin-top: -16px; }

span.preheader {
  display: none !important;
  visibility: hidden;
  mso-hide: all !important;
  font-size: 1px;
  color: #f3f3f3;
  line-height: 1px;
  max-height: 0px;
  max-width: 0px;
  opacity: 0;
  overflow: hidden; }

table.button {
  width: auto;
  margin: 0 0 16px 0;
  Margin: 0 0 16px 0; }
  table.button table td {
    text-align: left;
    color: #fefefe;
    background: #386DB0;
    border-radius: 8px;
    }
    table.button table td a {
      font-family: Helvetica, Arial, sans-serif;
      font-size: 16px;
      font-weight: lighter;
      color: #fefefe;
      text-decoration: none;
      display: inline-block;
      padding: 8px 16px 8px 16px;
     }
  table.button.radius table td {
    border-radius: 8px;
    border: none; }
  table.button.rounded table td {
    border-radius: 500px;
    border: none; }

table.button:hover table tr td a,
table.button:active table tr td a,
table.button table tr td a:visited,
table.button.tiny:hover table tr td a,
table.button.tiny:active table tr td a,
table.button.tiny table tr td a:visited,
table.button.small:hover table tr td a,
table.button.small:active table tr td a,
table.button.small table tr td a:visited,
table.button.large:hover table tr td a,
table.button.large:active table tr td a,
table.button.large table tr td a:visited {
  color: #fefefe; }

table.button.tiny table td,
table.button.tiny table a {
  padding: 4px 8px 4px 8px; }

table.button.tiny table a {
  font-size: 10px;
  font-weight: normal; }

table.button.small table td,
table.button.small table a {
  padding: 5px 10px 5px 10px;
  font-size: 12px; }

table.button.large table a {
  padding: 16px 32px 16px 32px;
  font-size: 20px; }

table.button.expand,
table.button.expanded {
  width: 100% !important; }
  table.button.expand table,
  table.button.expanded table {
    width: 100%; }
    table.button.expand table a,
    table.button.expanded table a {
      text-align: center;
      width: 100%;
      padding-left: 0;
      padding-right: 0; }
  table.button.expand center,
  table.button.expanded center {
    min-width: 0; }

table.button:hover table td,
table.button:visited table td,
table.button:active table td {
  background: #5689C9;
  color: #fefefe; }

table.button:hover table a,
table.button:visited table a,
table.button:active table a {
  border: 0 solid #5689C9; }

table.button.secondary table td {
  background: #777777;
  color: #fefefe;
  border: 0px solid #777777; }

table.button.secondary table a {
  color: #fefefe;
  border: 0 solid #777777; }

table.button.secondary:hover table td {
  background: #919191;
  color: #fefefe; }

table.button.secondary:hover table a {
  border: 0 solid #919191; }

table.button.secondary:hover table td a {
  color: #fefefe; }

table.button.secondary:active table td a {
  color: #fefefe; }

table.button.secondary table td a:visited {
  color: #fefefe; }

table.button.success table td {
  background: #3adb76;
  border: 0px solid #3adb76; }

table.button.success table a {
  border: 0 solid #3adb76; }

table.button.success:hover table td {
  background: #23bf5d; }

table.button.success:hover table a {
  border: 0 solid #23bf5d; }

table.button.alert table td {
  background: #ec5840;
  border: 0px solid #ec5840; }

table.button.alert table a {
  border: 0 solid #ec5840; }

table.button.alert:hover table td {
  background: #e23317; }

table.button.alert:hover table a {
  border: 0 solid #e23317; }

table.button.warning table td {
  background: #ffae00;
  border: 0px solid #ffae00; }

table.button.warning table a {
  border: 0px solid #ffae00; }

table.button.warning:hover table td {
  background: #cc8b00; }

table.button.warning:hover table a {
  border: 0px solid #cc8b00; }

table.callout {
  margin-bottom: 16px;
  Margin-bottom: 16px; }

th.callout-inner {
  width: 100%;
  border: 1px solid #cbcbcb;
  padding: 10px;
  background: #fefefe; }
  th.callout-inner.primary {
    background: #def0fc;
    border: 1px solid #444444;
    color: #0a0a0a; }
  th.callout-inner.secondary {
    background: #ebebeb;
    border: 1px solid #444444;
    color: #0a0a0a; }
  th.callout-inner.success {
    background: #e1faea;
    border: 1px solid #1b9448;
    color: #fefefe; }
  th.callout-inner.warning {
    background: #fff3d9;
    border: 1px solid #996800;
    color: #fefefe; }
  th.callout-inner.alert {
    background: #fce6e2;
    border: 1px solid #b42912;
    color: #fefefe; }

.thumbnail {
  border: solid 4px #fefefe;
  box-shadow: 0 0 0 1px rgba(10, 10, 10, 0.2);
  display: inline-block;
  line-height: 0;
  max-width: 100%;
  transition: box-shadow 200ms ease-out;
  border-radius: 3px;
  margin-bottom: 16px; }
  .thumbnail:hover, .thumbnail:focus {
    box-shadow: 0 0 6px 1px rgba(33, 153, 232, 0.5); }

table.menu {
  width: 580px; }
  table.menu td.menu-item,
  table.menu th.menu-item {
    padding: 10px;
    padding-right: 10px; }
    table.menu td.menu-item a,
    table.menu th.menu-item a {
      color: #2199e8; }

table.menu.vertical td.menu-item,
table.menu.vertical th.menu-item {
  padding: 10px;
  padding-right: 0;
  display: block; }
  table.menu.vertical td.menu-item a,
  table.menu.vertical th.menu-item a {
    width: 100%; }

table.menu.vertical td.menu-item table.menu.vertical td.menu-item,
table.menu.vertical td.menu-item table.menu.vertical th.menu-item,
table.menu.vertical th.menu-item table.menu.vertical td.menu-item,
table.menu.vertical th.menu-item table.menu.vertical th.menu-item {
  padding-left: 10px; }

table.menu.text-center a {
  text-align: center; }

.menu[align="center"] {
  width: auto !important; }

body.outlook p {
  display: inline !important; }

@media only screen and (max-width: 596px) {
  table.body img {
    width: auto;
    height: auto; }
  table.body center {
    min-width: 0 !important; }
  table.body .container {
    width: 95% !important; }
  table.body .columns,
  table.body .column {
    height: auto !important;
    -moz-box-sizing: border-box;
    -webkit-box-sizing: border-box;
    box-sizing: border-box;
    padding-left: 16px !important;
    padding-right: 16px !important; }
    table.body .columns .column,
    table.body .columns .columns,
    table.body .column .column,
    table.body .column .columns {
      padding-left: 0 !important;
      padding-right: 0 !important; }
  table.body .collapse .columns,
  table.body .collapse .column {
    padding-left: 0 !important;
    padding-right: 0 !important; }
  td.small-1,
  th.small-1 {
    display: inline-block !important;
    width: 8.33333% !important; }
  td.small-2,
  th.small-2 {
    display: inline-block !important;
    width: 16.66667% !important; }
  td.small-3,
  th.small-3 {
    display: inline-block !important;
    width: 25% !important; }
  td.small-4,
  th.small-4 {
    display: inline-block !important;
    width: 33.33333% !important; }
  td.small-5,
  th.small-5 {
    display: inline-block !important;
    width: 41.66667% !important; }
  td.small-6,
  th.small-6 {
    display: inline-block !important;
    width: 50% !important; }
  td.small-7,
  th.small-7 {
    display: inline-block !important;
    width: 58.33333% !important; }
  td.small-8,
  th.small-8 {
    display: inline-block !important;
    width: 66.66667% !important; }
  td.small-9,
  th.small-9 {
    display: inline-block !important;
    width: 75% !important; }
  td.small-10,
  th.small-10 {
    display: inline-block !important;
    width: 83.33333% !important; }
  td.small-11,
  th.small-11 {
    display: inline-block !important;
    width: 91.66667% !important; }
  td.small-12,
  th.small-12 {
    display: inline-block !important;
    width: 100% !important; }
  .columns td.small-12,
  .column td.small-12,
  .columns th.small-12,
  .column th.small-12 {
    display: block !important;
    width: 100% !important; }
  table.body td.small-offset-1,
  table.body th.small-offset-1 {
    margin-left: 8.33333% !important;
    Margin-left: 8.33333% !important; }
  table.body td.small-offset-2,
  table.body th.small-offset-2 {
    margin-left: 16.66667% !important;
    Margin-left: 16.66667% !important; }
  table.body td.small-offset-3,
  table.body th.small-offset-3 {
    margin-left: 25% !important;
    Margin-left: 25% !important; }
  table.body td.small-offset-4,
  table.body th.small-offset-4 {
    margin-left: 33.33333% !important;
    Margin-left: 33.33333% !important; }
  table.body td.small-offset-5,
  table.body th.small-offset-5 {
    margin-left: 41.66667% !important;
    Margin-left: 41.66667% !important; }
  table.body td.small-offset-6,
  table.body th.small-offset-6 {
    margin-left: 50% !important;
    Margin-left: 50% !important; }
  table.body td.small-offset-7,
  table.body th.small-offset-7 {
    margin-left: 58.33333% !important;
    Margin-left: 58.33333% !important; }
  table.body td.small-offset-8,
  table.body th.small-offset-8 {
    margin-left: 66.66667% !important;
    Margin-left: 66.66667% !important; }
  table.body td.small-offset-9,
  table.body th.small-offset-9 {
    margin-left: 75% !important;
    Margin-left: 75% !important; }
  table.body td.small-offset-10,
  table.body th.small-offset-10 {
    margin-left: 83.33333% !important;
    Margin-left: 83.33333% !important; }
  table.body td.small-offset-11,
  table.body th.small-offset-11 {
    margin-left: 91.66667% !important;
    Margin-left: 91.66667% !important; }
  table.body table.columns td.expander,
  table.body table.columns th.expander {
    display: none !important; }
  table.body .right-text-pad,
  table.body .text-pad-right {
    padding-left: 10px !important; }
  table.body .left-text-pad,
  table.body .text-pad-left {
    padding-right: 10px !important; }
  table.menu {
    width: 100% !important; }
    table.menu td,
    table.menu th {
      width: auto !important;
      display: inline-block !important; }
    table.menu.vertical td,
    table.menu.vertical th, table.menu.small-vertical td,
    table.menu.small-vertical th {
      display: block !important; }
  table.menu[align="center"] {
    width: auto !important; }
  table.button.small-expand,
  table.button.small-expanded {
    width: 100% !important; }
    table.button.small-expand table,
    table.button.small-expanded table {
      width: 100%; }
      table.button.small-expand table a,
      table.button.small-expanded table a {
        text-align: center !important;
        width: 100% !important;
        padding-left: 0 !important;
        padding-right: 0 !important; }
    table.button.small-expand center,
    table.button.small-expanded center {
      min-width: 0; } }

    </style>

    <style>
      body,
      html,
      .body {
        background: #f3f3f3 !important;
      }

      .header {
        background: #f3f3f3;
      }
    </style>
  </head>

  <body>
    <!-- <style> -->
    <table class="body" data-made-with-foundation="">
      <tr>
        <td class="float-center" align="center" valign="top">
          <center data-parsed="">
            <table class="spacer float-center">
              <tbody>
                <tr>
                  <td height="16px" style="font-size:16px;line-height:16px;"><br/></td>
                </tr>
              </tbody>
            </table>
            <table align="center" class="container float-center">
              <tbody>
                <tr>
                  <td>
                    <table class="row header">
                      <tbody>
                        <tr>
                          <th class="small-12 large-12 columns first last">
                            <table>
                              <tr>
                                <th>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="16px" style="font-size:16px;line-height:16px;"><br/></td>
                                      </tr>
                                    </tbody>
                                  </table>
                                  <h4 class="text-center"></h4> </th>
                                <th class="expander"></th>
                              </tr>
                            </table>
                          </th>
                        </tr>
                      </tbody>
                    </table>
                    <table class="row">
                      <tbody>
                        <tr>
                          <th class="small-12 large-12 columns first last">
                            <table>
                              <tr>
                                <th>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="32px" style="font-size:32px;line-height:32px;"><br/></td>
                                      </tr>
                                    </tbody>
                                  </table>
                                  <center data-parsed=""> <img src="../Logo.svg" align="center" class="float-left"> </center>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="16px" style="font-size:16px;line-height:16px;"><br/> </td>
                                      </tr>
                                    </tbody>
                                  </table>
                                  <h1 class="text-center">Hot off the presses!</h1>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="16px" style="font-size:16px;line-height:16px;"><br/></td>
                                      </tr>
                                    </tbody>
                                  </table>
                                  <p class="text-center">Your <a href="REPLACEME">dataset</a> is ready to be downloaded.</p>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="16px" style="font-size:16px;line-height:16px;"><br/></td>
                                      </tr>
                                    </tbody>
                                  </table>
                                  <table style="width:100%">
                                    <tr>
                                      <td align="center" valign="top">
                                        <table>
                                          <tr>
                                            <td align="center" valign="top" style="text-align:-webkit-center">
                                              <a href="REPLACEME#" align="center" class="float-center" style= " background-color:#386DB0; padding: 16px 32px 16px 32px;-webkit-border-radius: 3px; border-radius: 3px;color: #ffffff; display:inline-block;text-transform: uppercase;font-weight: lighter;">Download Now</a>
                                            </td>
                                          </tr>
                                        </table>
                                      </td>
                                      <td class="expander"></td>
                                    </tr>
                                  </table>
                                  <hr>
                                  <table class="spacer">
                                    <tbody>
                                      <tr>
                                        <td height="8px" style="font-size:8px;line-height:8px;"><br/></td>
                                      </tr>
                                    </tbody>
                                  </table>
                                    <p class="text-center"><small> <a href ="https://refine.bio/terms">Terms of Use</a></small></p>
                                    <p class="text-center"><small> Alex's Lemonade Stand Foundation <br/> 111 Presidential Blvd, Suite 203, Bala Cynwyd,  PA 19004</small></p>
                                </th>
                                <th class="expander"></th>
                              </tr>
                            </table>
                          </th>
                        </tr>
                      </tbody>
                    </table>
                    <table class="spacer">
                      <tbody>
                        <tr>
                          <td height="16px" style="font-size:16px;line-height:16px;"><br/></td>
                        </tr>
                      </tbody>
                    </table>
                  </td>
                </tr>
              </tbody>
            </table>
          </center>
        </td>
      </tr>
    </table>
  </body>

</html>
"""
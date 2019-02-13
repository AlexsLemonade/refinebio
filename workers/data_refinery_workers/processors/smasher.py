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

from botocore.exceptions import ClientError
from datetime import datetime, timedelta
from django.conf import settings
from django.utils import timezone
from pathlib import Path
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr
from sklearn import preprocessing
from typing import Dict
import numpy as np
import pandas as pd

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
logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """
    Fetches and prepares the files to smash.
    """

    all_sample_files = []
    job_context['input_files'] = {}

    # `key` can either be the species name or experiment accession.
    for key, samples in job_context["samples"].items():
        smashable_files = []
        for sample in samples:
            smashable_file = sample.get_most_recent_smashable_result_file()

            if smashable_file is not None:
                smashable_files = smashable_files + [smashable_file]
        smashable_files = list(set(smashable_files))
        job_context['input_files'][key] = smashable_files
        all_sample_files = all_sample_files + smashable_files

    # Filter empty results. This shouldn't get here, but it's possible, so we filter just in case it does.
    all_sample_files = [sf for sf in all_sample_files if sf is not None]

    if all_sample_files == []:
        error_message = "Couldn't get any files to smash for Smash job!!"
        logger.error(error_message,
                     dataset_id=job_context['dataset'].id,
                     samples=job_context["samples"])

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

    return job_context


def _add_annotation_column(annotation_columns, column_name):
    """Add annotation column names in place.
    Any column_name that starts with "refinebio_" will be skipped.
    """

    if not column_name.startswith("refinebio_"):
        annotation_columns.add(column_name)


def _get_tsv_columns(samples_metadata):
    """Returns an array of strings that will be written as a TSV file's
    header. The columns are based on fields found in samples_metadata.

    Some nested annotation fields are taken out as separate columns
    because they are more important than the others.
    """

    refinebio_columns = set()
    annotation_columns = set()
    for sample_metadata in samples_metadata.values():
        for meta_key, meta_value in sample_metadata.items():
            if meta_key != 'refinebio_annotations':
                refinebio_columns.add(meta_key)
                continue

            # Decompose sample_metadata["annotations"], which is an array of annotations!
            for annotation in meta_value:
                for annotation_key, annotation_value in annotation.items():
                    # For ArrayExpress samples, take out the fields
                    # nested in "characteristic" as separate columns.
                    if (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                        and annotation_key == "characteristic"):
                        for pair_dict in annotation_value:
                            if 'category' in pair_dict and 'value' in pair_dict:
                                _add_annotation_column(annotation_columns, pair_dict['category'])
                    # For ArrayExpress samples, also take out the fields
                    # nested in "variable" as separate columns.
                    elif (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                          and annotation_key == "variable"):
                        for pair_dict in annotation_value:
                            if 'name' in pair_dict and 'value' in pair_dict:
                                _add_annotation_column(annotation_columns, pair_dict['name'])
                    # For ArrayExpress samples, skip "source" field
                    elif (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                          and annotation_key == "source"):
                        continue
                    # For GEO samples, take out the fields nested in
                    # "characteristics_ch1" as separate columns.
                    elif (sample_metadata.get('refinebio_source_database', '') == "GEO"
                          and annotation_key == "characteristics_ch1"): # array of strings
                        for pair_str in annotation_value:
                            if ':' in pair_str:
                                tokens = pair_str.split(':', 1)
                                _add_annotation_column(annotation_columns, tokens[0])
                    # Saves all other annotation fields in separate columns
                    else:
                        _add_annotation_column(annotation_columns, annotation_key)

    # Return sorted columns, in which "refinebio_accession_code" is always the first,
    # followed by the other refinebio columns (in alphabetic order), and
    # annotation columns (in alphabetic order) at the end.
    refinebio_columns.discard('refinebio_accession_code')
    return ['refinebio_accession_code'] + sorted(refinebio_columns) + sorted(annotation_columns)

def _add_annotation_value(row_data, col_name, col_value, sample_accession_code):
    """Adds a new `col_name` key whose value is `col_value` to row_data.
    If col_name already exists in row_data with different value, print
    out a warning message.
    """

    # Generate a warning message if annotation field name starts with
    # "refinebio_".  This should rarely (if ever) happen.
    if col_name.startswith("refinebio_"):
        logger.warning(
            "Annotation value skipped",
            annotation_field=col_name,
            annotation_value=col_value,
            sample_accession_code=sample_accession_code
        )
    elif col_name not in row_data:
        row_data[col_name] = col_value
    # Generate a warning message in case of conflicts of annotation values.
    # (Requested by Dr. Jackie Taroni)
    elif row_data[col_name] != col_value:
        logger.warning(
            "Conflict of values found in column %s: %s vs. %s" % (
                col_name, row_data[col_name], col_value),
            sample_accession_code=sample_accession_code
        )


def _get_tsv_row_data(sample_metadata):
    """Returns field values based on input sample_metadata.

    Some annotation fields are treated specially because they are more
    important.  See `_get_tsv_columns` function above for details.
    """

    sample_accession_code = sample_metadata.get('refinebio_accession_code', '')
    row_data = dict()
    for meta_key, meta_value in sample_metadata.items():
        # If the field is a refinebio-specific field, simply copy it.
        if meta_key != 'refinebio_annotations':
            row_data[meta_key] = meta_value
            continue

        # Decompose sample_metadata["refinebio_annotations"], which is
        # an array of annotations.
        for annotation in meta_value:
            for annotation_key, annotation_value in annotation.items():
                # "characteristic" in ArrayExpress annotation
                if (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                    and annotation_key == "characteristic"):
                    for pair_dict in annotation_value:
                        if 'category' in pair_dict and 'value' in pair_dict:
                            col_name, col_value = pair_dict['category'], pair_dict['value']
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                # "variable" in ArrayExpress annotation
                elif (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                      and annotation_key == "variable"):
                    for pair_dict in annotation_value:
                        if 'name' in pair_dict and 'value' in pair_dict:
                            col_name, col_value = pair_dict['name'], pair_dict['value']
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                 # Skip "source" field ArrayExpress sample's annotation
                elif (sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS"
                      and annotation_key == "source"):
                    continue
                # "characteristics_ch1" in GEO annotation
                elif (sample_metadata.get('refinebio_source_database', '') == "GEO"
                      and annotation_key == "characteristics_ch1"): # array of strings
                    for pair_str in annotation_value:
                        if ':' in pair_str:
                            col_name, col_value = pair_str.split(':', 1)
                            col_value = col_value.strip()
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                # If annotation_value includes only a 'name' key, extract its value directly:
                elif (isinstance(annotation_value, dict)
                      and len(annotation_value) == 1 and 'name' in annotation_value):
                    _add_annotation_value(row_data, annotation_key, annotation_value['name'],
                                          sample_accession_code)
                # If annotation_value is a single-element array, extract the element directly:
                elif isinstance(annotation_value, list) and len(annotation_value) == 1:
                    _add_annotation_value(row_data, annotation_key, annotation_value[0],
                                          sample_accession_code)
                # Otherwise save all annotation fields in separate columns
                else:
                    _add_annotation_value(row_data, annotation_key, annotation_value,
                                          sample_accession_code)

    return row_data


def _write_tsv_json(job_context, metadata, smash_path):
    """Writes tsv files on disk.
    If the dataset is aggregated by species, also write species-level
    JSON file.
    """

    # Uniform TSV header per dataset
    columns = _get_tsv_columns(metadata['samples'])

    # Per-Experiment Metadata
    if job_context["dataset"].aggregate_by == "EXPERIMENT":
        tsv_paths = []
        for experiment_title, experiment_data in metadata['experiments'].items():
            experiment_dir = smash_path + experiment_title + '/'
            experiment_dir = experiment_dir.encode('ascii', 'ignore')
            os.makedirs(experiment_dir, exist_ok=True)
            tsv_path = experiment_dir.decode("utf-8") + 'metadata_' + experiment_title + '.tsv'
            tsv_path = tsv_path.encode('ascii', 'ignore')
            tsv_paths.append(tsv_path)
            with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
                dw.writeheader()
                for sample_title, sample_metadata in metadata['samples'].items():
                    if sample_title in experiment_data['sample_titles']:
                        row_data = _get_tsv_row_data(sample_metadata)
                        dw.writerow(row_data)
        return tsv_paths
    # Per-Species Metadata
    elif job_context["dataset"].aggregate_by == "SPECIES":
        tsv_paths = []
        for species in job_context['input_files'].keys():
            species_dir = smash_path + species + '/'
            os.makedirs(species_dir, exist_ok=True)
            samples_in_species = []
            tsv_path = species_dir + "metadata_" + species + '.tsv'
            tsv_paths.append(tsv_path)
            with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
                dw.writeheader()
                for sample_metadata in metadata['samples'].values():
                    if sample_metadata.get('refinebio_organism', '') == species:
                        row_data = _get_tsv_row_data(sample_metadata)
                        dw.writerow(row_data)
                        samples_in_species.append(sample_metadata)

            # Writes a json file for current species:
            if len(samples_in_species):
                species_metadata = {
                    'species': species,
                    'samples': samples_in_species
                }
                json_path = species_dir + "metadata_" + species + '.json'
                with open(json_path, 'w', encoding='utf-8') as json_file:
                    json.dump(species_metadata, json_file, indent=4, sort_keys=True)
        return tsv_paths
    # All Metadata
    else:
        all_dir = smash_path + "ALL/"
        os.makedirs(all_dir, exist_ok=True)
        tsv_path = all_dir + 'metadata_ALL.tsv' 
        with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
            dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
            dw.writeheader()
            for sample_metadata in metadata['samples'].values():
                row_data = _get_tsv_row_data(sample_metadata)
                dw.writerow(row_data)
        return [tsv_path]

def _quantile_normalize(job_context: Dict, ks_check=True, ks_stat=0.001) -> Dict:
    """
    Apply quantile normalization.

    """
    # Prepare our QN target file
    organism = job_context['organism']
    qn_target = utils.get_most_recent_qn_target_for_organism(organism)

    if not qn_target:
        logger.error("Could not find QN target for Organism!",
            organism=organism,
            dataset_id=job_context['dataset'].id,
            dataset_data=job_context['dataset'].data,
            processor_job_id=job_context["job"].id,
        )
        job_context['dataset'].success = False
        job_context['job'].failure_reason = "Could not find QN target for Organism: " + str(organism)
        job_context['dataset'].failure_reason = "Could not find QN target for Organism: " + str(organism)
        job_context['dataset'].save()
        job_context['job'].success = False
        job_context['failure_reason'] = "Could not find QN target for Organism: " + str(organism)
        return job_context
    else:
        qn_target_path = qn_target.sync_from_s3()
        qn_target_frame = pd.read_csv(qn_target_path, sep='\t', header=None,
                                      index_col=None, error_bad_lines=False)

        # Prepare our RPy2 bridge
        pandas2ri.activate()
        preprocessCore = importr('preprocessCore')
        as_numeric = rlang("as.numeric")
        data_matrix = rlang('data.matrix')

        # Convert the smashed frames to an R numeric Matrix
        # and the target Dataframe into an R numeric Vector
        target_vector = as_numeric(qn_target_frame[0])
        merged_matrix = data_matrix(job_context['merged_no_qn'])

        # Perform the Actual QN
        reso = preprocessCore.normalize_quantiles_use_target(
                                            x=merged_matrix,
                                            target=target_vector,
                                            copy=True
                                        )

        # Verify this QN, related: https://github.com/AlexsLemonade/refinebio/issues/599#issuecomment-422132009
        set_seed = rlang("set.seed")
        combn = rlang("combn")
        ncol = rlang("ncol")
        ks_test = rlang("ks.test")
        which = rlang("which")

        set_seed(123)

        n = ncol(reso)[0]
        m = 2
        if n < m:
            raise Exception("Found fewer columns than required for QN combinatorial - bad smash?")
        combos = combn(ncol(reso), 2)

        # Convert to NP, Shuffle, Return to R
        ar = np.array(combos)
        np.random.shuffle(np.transpose(ar))
        nr, nc = ar.shape
        combos = ro.r.matrix(ar, nrow=nr, ncol=nc)

        # adapted from
        # https://stackoverflow.com/questions/9661469/r-t-test-over-all-columns
        # apply KS test to randomly selected pairs of columns (samples)
        for i in range(1, min(ncol(combos)[0], 100)):
            value1 = combos.rx(1, i)[0]
            value2 = combos.rx(2, i)[0]

            test_a = reso.rx(True, value1)
            test_b = reso.rx(True, value2)

            # RNA-seq has a lot of zeroes in it, which
            # breaks the ks_test. Therefore we want to
            # filter them out. To do this we drop the
            # lowest half of the values. If there's
            # still zeroes in there, then that's
            # probably too many zeroes so it's okay to
            # fail.
            median_a = np.median(test_a)
            median_b = np.median(test_b)

            # `which` returns indices which are
            # 1-indexed. Python accesses lists with
            # zero-indexes, even if that list is
            # actually an R vector. Therefore subtract
            # 1 to account for the difference.
            test_a = [test_a[i-1] for i in which(test_a > median_a)]
            test_b = [test_b[i-1] for i in which(test_b > median_b)]

            # The python list comprehension gives us a
            # python list, but ks_test wants an R
            # vector so let's go back.
            test_a = as_numeric(test_a)
            test_b = as_numeric(test_b)

            ks_res = ks_test(test_a, test_b)
            statistic = ks_res.rx('statistic')[0][0]
            pvalue = ks_res.rx('p.value')[0][0]

            job_context['ks_statistic'] = statistic
            job_context['ks_pvalue'] = pvalue

            # We're unsure of how strigent to be about
            # the pvalue just yet, so we're extra lax
            # rather than failing tons of tests. This may need tuning.
            if ks_check:
                if statistic > ks_stat or pvalue < 0.8:
                    raise Exception("Failed Kolmogorov Smirnov test! Stat: " +
                                    str(statistic) + ", PVal: " + str(pvalue))

        # And finally convert back to Pandas
        ar = np.array(reso)
        new_merged = pd.DataFrame(ar, columns=job_context['merged_no_qn'].columns, index=job_context['merged_no_qn'].index)
        job_context['merged_qn'] = new_merged
        merged = new_merged
    return job_context   

def _smash(job_context: Dict, how="inner") -> Dict:
    """
    Smash all of the samples together!

    Steps:
        Combine common genes (pandas merge)
        Transpose such that genes are columns (features)
        Scale features with sci-kit learn
        Transpose again such that samples are columns and genes are rows
    """

    # We have already failed - return now so we can send our fail email.
    if job_context['dataset'].failure_reason not in ['', None]:
        return job_context

    try:
        # Prepare the output directory
        smash_path = job_context["output_dir"]

        scalers = {
            'MINMAX': preprocessing.MinMaxScaler,
            'STANDARD': preprocessing.StandardScaler,
            'ROBUST': preprocessing.RobustScaler,
        }

        unsmashable_files = []
        num_samples = 0

        # Smash all of the sample sets
        logger.debug("About to smash!",
                     input_files=job_context['input_files'],
                     dataset_data=job_context['dataset'].data,
        )

        job_context['technologies'] = {'microarray': [], 'rnaseq': []}
        # Once again, `key` is either a species name or an experiment accession
        for key, input_files in job_context['input_files'].items():

            # Merge all the frames into one
            all_frames = []

            for computed_file in input_files:

                try:
                    # Download the file to a job-specific location so it
                    # won't disappear while we're using it.

                    computed_file_path = job_context["work_dir"] + computed_file.filename
                    computed_file_path = computed_file.get_synced_file_path(path=computed_file_path)

                    # Bail appropriately if this isn't a real file.
                    if not computed_file_path or not os.path.exists(computed_file_path):
                        unsmashable_files.append(computed_file_path)
                        logger.error("Smasher received non-existent file path.",
                            computed_file_path=computed_file_path,
                            computed_file=computed_file,
                            dataset=job_context['dataset'],
                            )
                        continue

                    data = _load_and_sanitize_file(computed_file_path)

                    if len(data.columns) > 2:
                        # Most of the time, >1 is actually bad, but we also need to support
                        # two-channel samples. I think ultimately those should be given some kind of
                        # special consideration.
                        logger.info("Found a frame with more than 2 columns - this shouldn't happen!",
                            computed_file_path=computed_file_path,
                            computed_file_id=computed_file.id
                            )
                        continue

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
                        data.columns = [computed_file.samples.all()[0].accession_code]
                    except ValueError:
                        # This sample might have multiple channels, or something else.
                        # Don't mess with it.
                        pass
                    except Exception as e:
                        # Okay, somebody probably forgot to create a SampleComputedFileAssociation
                        data.columns = [computed_file.filename]

                    if computed_file_path.endswith("lengthScaledTPM.tsv"):
                        job_context['technologies']['rnaseq'].append(data.columns)
                    else:
                        job_context['technologies']['microarray'].append(data.columns)

                    all_frames.append(data)
                    num_samples = num_samples + 1

                    if (num_samples % 100) == 0:
                        logger.info("Loaded " + str(num_samples) + " samples into frames.",
                            dataset_id=job_context['dataset'].id,
                            how=how
                        )

                except Exception as e:
                    unsmashable_files.append(computed_file_path)
                    logger.exception("Unable to smash file",
                        file=computed_file_path,
                        dataset_id=job_context['dataset'].id,
                        )
                finally:
                    # Delete before archiving the work dir
                    if computed_file_path:
                        os.remove(computed_file_path)

            job_context['all_frames'] = all_frames

            if len(all_frames) < 1:
                logger.warning("Was told to smash a frame with no frames!",
                    key=key,
                    input_files=str(input_files)
                )
                continue

            # Merge all of the frames we've gathered into a single big frame, skipping duplicates.
            # TODO: If the very first frame is the wrong platform, are we boned?
            merged = all_frames[0]
            i = 1

            old_len_merged = len(merged)
            new_len_merged = len(merged)
            merged_backup = merged

            if how == "inner":
                while i < len(all_frames):
                    frame = all_frames[i]
                    i = i + 1

                    if i % 1000 == 0:
                        logger.info("Smashing keyframe",
                            i=i
                        )

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
                                       processor_job_id=job_context["job"].id,
                                       column=column,
                                       frame=frame
                        )
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
                            old_len_merged=old_len_merged,
                            new_len_merged=new_len_merged,
                            bad_frame_number=i,
                        )
                        merged = merged_backup
                        new_len_merged = len(merged)
                        try:
                            unsmashable_files.append(frame.columns[0])
                        except Exception:
                            # Something is really, really wrong with this frame.
                            pass

                    old_len_merged = len(merged)
                    merged_backup = merged
            else:
                merged = pd.concat(all_frames, axis=1, keys=None, join='outer', copy=False)

            job_context['original_merged'] = merged

            # Quantile Normalization
            if job_context['dataset'].quantile_normalize:
                try:
                    job_context['merged_no_qn'] = merged
                    job_context['organism'] = computed_file.samples.first().organism
                    job_context = _quantile_normalize(job_context)
                    merged = job_context.get('merged_qn', None)
                    # We probably don't have an QN target or there is another error,
                    # so let's fail gracefully.
                    if merged is None:
                        e = "Problem occured during quantile normalization: No merged_qn"
                        logger.error(e,
                            dataset_id=job_context['dataset'].id,
                            dataset_data=job_context['dataset'].data,
                            processor_job_id=job_context["job"].id,
                        )
                        job_context['dataset'].success = False
                        job_context['job'].failure_reason = "Failure reason: " + str(e)
                        job_context['dataset'].failure_reason = "Failure reason: " + str(e)
                        job_context['dataset'].save()
                        # Delay failing this pipeline until the failure notify has been sent
                        job_context['job'].success = False
                        job_context['failure_reason'] = str(e)
                        return job_context
                except Exception as e:
                    logger.exception("Problem occured during quantile normalization",
                        dataset_id=job_context['dataset'].id,
                        dataset_data=job_context['dataset'].data,
                        processor_job_id=job_context["job"].id,
                    )
                    job_context['dataset'].success = False
                    job_context['job'].failure_reason = "Failure reason: " + str(e)
                    job_context['dataset'].failure_reason = "Failure reason: " + str(e)
                    job_context['dataset'].save()
                    # Delay failing this pipeline until the failure notify has been sent
                    job_context['job'].success = False
                    job_context['failure_reason'] = str(e)
                    return job_context
            # End QN

            # Transpose before scaling
            # Do this even if we don't want to scale in case transpose
            # modifies the data in any way. (Which it shouldn't but
            # we're paranoid.)
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
            subdir = ''
            if job_context['dataset'].aggregate_by in ["SPECIES", "EXPERIMENT"]:
                subdir = key
            elif job_context['dataset'].aggregate_by == "ALL":
                subdir = "ALL"

            outfile_dir = smash_path + key + "/"
            os.makedirs(outfile_dir, exist_ok=True)
            outfile = outfile_dir + key + ".tsv"
            job_context['smash_outfile'] = outfile
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

        # Write samples metadata to TSV
        try:
            tsv_paths = _write_tsv_json(job_context, metadata, smash_path)
            job_context['metadata_tsv_paths'] = tsv_paths
            # Metadata to JSON
            metadata['created_at'] = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')
            with open(smash_path + 'aggregated_metadata.json', 'w', encoding='utf-8') as metadata_file:
                json.dump(metadata, metadata_file, indent=4, sort_keys=True)
        except Exception as e:
            logger.exception("Failed to write metadata TSV!",
                j_id = job_context['job'].id)
            job_context['metadata_tsv_paths'] = None
        metadata['files'] = os.listdir(smash_path)

        # Finally, compress all files into a zip
        final_zip_base = "/home/user/data_store/smashed/" + str(job_context["dataset"].pk)
        shutil.make_archive(final_zip_base, 'zip', smash_path)
        job_context["output_file"] = final_zip_base + ".zip"
    except Exception as e:
        logger.exception("Could not smash dataset.",
                        dataset_id=job_context['dataset'].id,
                        processor_job_id=job_context['job_id'],
                        input_files=job_context['input_files'])
        job_context['dataset'].success = False
        job_context['job'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].failure_reason = "Failure reason: " + str(e)
        job_context['dataset'].save()
        # Delay failing this pipeline until the failure notify has been sent
        job_context['job'].success = False
        job_context['failure_reason'] = str(e)
        return job_context

    job_context['metadata'] = metadata
    job_context['unsmashable_files'] = unsmashable_files
    job_context['dataset'].success = True
    job_context['dataset'].save()

    logger.debug("Created smash output!",
        archive_location=job_context["output_file"])

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

def _upload(job_context: Dict) -> Dict:
    """ Uploads the result file to S3 and notifies user. """

    # There has been a failure already, don't try to upload anything.
    if not job_context.get("output_file", None):
        logger.error("Was told to upload a smash result without an output_file.")
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

        # Don't send an email if we don't have address.
        if job_context["dataset"].email_address:
            SENDER = "Refine.bio Mail Robot <noreply@refine.bio>"
            RECIPIENT = job_context["dataset"].email_address
            AWS_REGION = "us-east-1"
            CHARSET = "UTF-8"

            # Link to the dataset page, where the user can re-try the download job
            dataset_url = 'https://www.refine.bio/dataset/' + str(job_context['dataset'].id)

            if job_context['job'].failure_reason not in ['', None]:
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
                BODY_TEXT = "Hot off the presses:\n\n" + job_context["result_url"] + "\n\nLove!,\nThe refine.bio Team"
                FORMATTED_HTML = BODY_HTML.replace('REPLACE_DOWNLOAD_URL', job_context["result_url"])\
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
    """Closes out the dataset object."""

    dataset = job_context["dataset"]
    dataset.is_processing = False
    dataset.is_processed = True
    dataset.is_available = True
    dataset.expires_on = timezone.now() + timedelta(days=1)
    dataset.save()

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
                        utils.end_job])

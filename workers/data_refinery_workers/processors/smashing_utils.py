# -*- coding: utf-8 -*-

import csv
import logging
import math
import multiprocessing
import os
import shutil
import time
from pathlib import Path
from typing import Dict, List, Tuple

from django.utils import timezone
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr
import numpy as np
import pandas as pd
import psutil
import rpy2.robjects as ro
import simplejson as json

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputedFile, Sample
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

# Take one fewer than 1/2 the total available threads
# also make the minimum threads 1.
# Use floor here because multiprocessing raises an exception if this isn't an int.
MULTIPROCESSING_WORKER_COUNT = max(1, math.floor(multiprocessing.cpu_count()/2) - 1)
MULTIPROCESSING_CHUNK_SIZE = 2000
RESULTS_BUCKET = get_env_variable("S3_RESULTS_BUCKET_NAME", "refinebio-results-bucket")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
BODY_HTML = Path(
    'data_refinery_workers/processors/smasher_email.min.html'
).read_text().replace('\n', '')
BODY_ERROR_HTML = Path(
    'data_refinery_workers/processors/smasher_email_error.min.html'
).read_text().replace('\n', '')
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


def log_state(message, job_id, start_time=False):
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process(os.getpid())
        ram_in_GB = process.memory_info().rss / BYTES_IN_GB
        logger.debug(message,
                     total_cpu=psutil.cpu_percent(),
                     process_ram=ram_in_GB,
                     job_id=job_id)

        if start_time:
            logger.debug('Duration: %s' % (time.time() - start_time), job_id=job_id)
        else:
            return time.time()


def prepare_files(job_context: Dict) -> Dict:
    """
    Fetches and prepares the files to smash.
    """
    start_prepare_files = log_state("start prepare files", job_context["job"].id)
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

    job_context['num_input_files'] = len(job_context['input_files'])
    job_context['group_by_keys'] = list(job_context['input_files'].keys())

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
        job_context["job"].failure_reason = ("Couldn't get any files to smash for Smash job"
                                             " - empty all_sample_files")
        return job_context

    dataset_id = str(job_context["dataset"].pk)
    job_context["work_dir"] = "/home/user/data_store/smashed/" + dataset_id + "/"
    # Ensure we have a fresh smash directory
    shutil.rmtree(job_context["work_dir"], ignore_errors=True)
    os.makedirs(job_context["work_dir"])

    job_context["output_dir"] = job_context["work_dir"] + "output/"
    os.makedirs(job_context["output_dir"])
    log_state("end prepare files", job_context["job"].id, start_prepare_files)
    return job_context


def _load_and_sanitize_file(computed_file_path, index=None):
    """ Read and sanitize a computed file """

    data = pd.read_csv(computed_file_path,
                       sep='\t',
                       header=0,
                       index_col=0,
                       error_bad_lines=False)

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
    #       fgenesh2_kg.7__3016__AT5G35080.1 (via http://plants.ensembl.org/Arabidopsis_lyrata/ \
    #       Gene/Summary?g=fgenesh2_kg.7__3016__AT5G35080.1;r=7:17949732-17952000;t=fgenesh2_kg. \
    #       7__3016__AT5G35080.1;db=core)
    data.index = data.index.str.replace(r"(\.[^.]*)$", '')

    # Squish duplicated rows together.
    # XXX/TODO: Is mean the appropriate method here?
    #           We can make this an option in future.
    # Discussion here: https://github.com/AlexsLemonade/refinebio/issues/186#issuecomment-395516419
    data = data.groupby(data.index, sort=False).mean()

    if index is not None:
        data = data.reindex(index)

    return data


def process_frame(work_dir,
                  computed_file,
                  sample_accession_code,
                  dataset_id,
                  aggregate_by,
                  index,
                  gene_ids,
                  job_id) -> Dict:

    log_state('processing frame {}'.format(index), job_id)

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
        computed_file_path = computed_file.get_synced_file_path(
            path="%s%s" % (work_dir, computed_file.filename)
        )

        # Bail appropriately if this isn't a real file.
        if not computed_file_path or not os.path.exists(computed_file_path):
            logger.warning("Smasher received non-existent file path.",
                           computed_file_path=computed_file_path,
                           computed_file_id=computed_file.id,
                           dataset_id=dataset_id)
            return unsmashable(computed_file.filename)

        data = _load_and_sanitize_file(computed_file_path, gene_ids)

        if len(data.columns) > 2:
            # Most of the time, >1 is actually bad, but we also need to support
            # two-channel samples. I think ultimately those should be given some kind of
            # special consideration.
            logger.info("Found a frame with more than 2 columns - this shouldn't happen!",
                        computed_file_path=computed_file_path,
                        computed_file_id=computed_file.id)
            return unsmashable(computed_file_path)

        # via https://github.com/AlexsLemonade/refinebio/issues/330:
        #   aggregating by experiment -> return untransformed output from tximport
        #   aggregating by species -> log2(x + 1) tximport output
        if aggregate_by == 'SPECIES' \
           and computed_file_path.endswith("lengthScaledTPM.tsv"):
            data = data + 1
            data = np.log2(data)

        # Detect if this data hasn't been log2 scaled yet.
        # Ideally done in the NO-OPPER, but sanity check here.
        if (not computed_file_path.endswith("lengthScaledTPM.tsv")) and (data.max() > 100).any():
            logger.info("Detected non-log2 microarray data.", computed_file_id=computed_file.id)
            data = np.log2(data)

        # Explicitly title this dataframe
        try:
            # Unfortuantely, we can't use this as `title` can cause a collision
            # data.columns = [computed_file.samples.all()[0].title]
            # So we use this, which also helps us support the case of missing
            # SampleComputedFileAssociation
            data.columns = [sample_accession_code]
        except ValueError as e:
            # This sample might have multiple channels, or something else.
            # Don't mess with it.
            logger.warn("Smasher found multi-channel column (probably) - skipping!",
                        exc_info=1,
                        computed_file_path=computed_file_path,)
            return unsmashable(computed_file.filename)
        except Exception as e:
            # Okay, somebody probably forgot to create a SampleComputedFileAssociation
            # Don't mess with it.
            logger.warn("Smasher found very bad column title - skipping!",
                        exc_info=1,
                        computed_file_path=computed_file_path)
            return unsmashable(computed_file.filename)

        is_rnaseq = computed_file_path.endswith("lengthScaledTPM.tsv")
        frame['technology'] = 'rnaseq' if is_rnaseq else 'microarray'

    except Exception as e:
        logger.exception("Unable to smash file",
                         file=computed_file_path,
                         dataset_id=dataset_id)
        return unsmashable(computed_file.filename)
    # TEMPORARY for iterating on compendia more quickly.
    # finally:
    #     # Delete before archiving the work dir
    #     if computed_file_path and os.path.exists(computed_file_path):
    #         os.remove(computed_file_path)

    return smashable(data)


def process_frames_for_key(key: str,
                           input_files: List[ComputedFile],
                           job_context: Dict,
                           merge_strategy: str) -> Dict:
    """Download, read, and chunk processed sample files from s3.

    `key` is the species or experiment whose samples are contained in `input_files`.

    Will populate add to job_context the keys 'microarray_frames' and
    'rnaseq_frames' with pandas dataframes containing
    MULTIPROCESSING_CHUNK_SIZE samples worth of data. Also adds the
    key 'unsmashable_files' containing a list of paths that were
    determined to be unsmashable.

    `merge_strategy` dictates how the chunks will be merged and must be
    one of the two values `inner` or `outer.
    """
    if merge_strategy != 'inner' and merge_strategy != 'outer':
        raise ValueError("merge_strategy must be either of the values 'inner' or 'outer'.")

    job_context['original_merged'] = pd.DataFrame()

    start_gene_ids = log_state("Collecting all gene identifiers for key {}".format(key),
                               job_context["job"].id)

    # Build up a list of gene identifiers because these will be the
    # rows of our matrices, and we want to preallocate them so we need
    # to know them all.
    all_gene_identifiers = set()

    microarray_columns = []
    rnaseq_columns = []
    for index, (computed_file, sample) in enumerate(input_files):
        frame = process_frame(job_context["work_dir"],
                              computed_file,
                              sample.accession_code,
                              job_context['dataset'].id,
                              job_context['dataset'].aggregate_by,
                              index,
                              None,
                              job_context["job"].id)

        # Count how many frames are in each tech so we can preallocate
        # the matrices in both directions.
        if not frame['unsmashable']:
            all_gene_identifiers = all_gene_identifiers.union(frame['dataframe'].index)

            # Each dataframe should only have 1 column, but it's returned as a list so use extend.
            if frame['technology'] == 'microarray':
                microarray_columns.extend(frame['dataframe'].columns)
            elif frame['technology'] == 'rnaseq':
                rnaseq_columns.extend(frame['dataframe'].columns)

    all_gene_identifiers = list(all_gene_identifiers)
    all_gene_identifiers.sort()

    log_template = ("Collected {0} gene identifiers for {1} across"
                    " {2} micrarry samples and {3} RNA-Seq samples.")
    log_state(log_template.format(len(all_gene_identifiers),
                                  key,
                                  len(microarray_columns),
                                  len(rnaseq_columns)),
              job_context["job"].id,
              start_gene_ids)
    start_build_matrix = log_state("Beginning to build the full matrices.",
                                   job_context["job"].id)

    # Preallocate the matrices to be the exact size we will need. This
    # should prevent any operations from happening while we build it
    # up, so the only RAM used will be needed.
    job_context['microarray_matrix'] = pd.DataFrame(data=None,
                                                    index=all_gene_identifiers,
                                                    columns=microarray_columns,
                                                    dtype=np.float64)
    job_context['rnaseq_matrix'] = pd.DataFrame(data=None,
                                                index=all_gene_identifiers,
                                                columns=rnaseq_columns,
                                                dtype=np.float64)

    for index, (computed_file, sample) in enumerate(input_files):
        processed_frame = process_frame(job_context["work_dir"],
                                        computed_file,
                                        sample.accession_code,
                                        job_context['dataset'].id,
                                        job_context['dataset'].aggregate_by,
                                        index,
                                        all_gene_identifiers,
                                        job_context["job"].id)

        if frame['unsmashable']:
            job_context['unsmashable_files'].append(frame['unsmashable_file'])
        else:
            # The dataframe for each sample will only have one column
            # whose header will be the accession code.
            column = processed_frame['dataframe'].columns[0]
            if processed_frame['technology'] == 'microarray':
                job_context['microarray_matrix'][column] = processed_frame['dataframe'].values
            elif processed_frame['technology'] == 'rnaseq':
                job_context['rnaseq_matrix'][column] = processed_frame['dataframe'].values

        del processed_frame

    job_context['num_samples'] = 0
    if job_context['microarray_matrix'] is not None:
        job_context['num_samples'] += len(job_context['microarray_matrix'].columns)
    if job_context['rnaseq_matrix'] is not None:
        job_context['num_samples'] += len(job_context['rnaseq_matrix'].columns)

    log_state("Built full matrices for key {}".format(key),
              job_context["job"].id,
              start_build_matrix)

    return job_context


def quantile_normalize(job_context: Dict, ks_check=True, ks_stat=0.001) -> Dict:
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
                     processor_job_id=job_context["job"].id,)
        job_context['dataset'].success = False
        failure_reason = "Could not find QN target for Organism: " + str(organism)
        job_context['job'].failure_reason = failure_reason
        job_context['dataset'].failure_reason = failure_reason
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

        # Verify this QN, related:
        # https://github.com/AlexsLemonade/refinebio/issues/599#issuecomment-422132009
        set_seed = rlang("set.seed")
        combn = rlang("combn")
        ncol = rlang("ncol")
        ks_test = rlang("ks.test")
        which = rlang("which")

        set_seed(123)

        n = ncol(reso)[0]
        m = 2
        if n >= m:
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
                if ks_check and (statistic > ks_stat or pvalue < 0.8):
                    job_context['ks_warning'] = ("Failed Kolmogorov Smirnov test! Stat: " +
                                                 str(statistic) + ", PVal: " + str(pvalue))
        else:
            logger.warning(("Not enough columns to perform KS test -"
                            " either bad smash or single saple smash."),
                           dataset_id=job_context['dataset'].id)

        # And finally convert back to Pandas
        ar = np.array(reso)
        new_merged = pd.DataFrame(ar,
                                  columns=job_context['merged_no_qn'].columns,
                                  index=job_context['merged_no_qn'].index)

        # Remove un-quantiled normalized matrix from job_context
        # because we no longer need it.
        job_context.pop('merged_no_qn')

        # And add the quantile normalized matrix to job_context.
        job_context['merged_qn'] = new_merged
    return job_context


def compile_metadata(job_context: Dict) -> Dict:
    """Compiles metadata about the job.

    Returns a new dict containing the metadata, not the job_context.
    """
    metadata = {}

    metadata['num_samples'] = job_context['num_samples']
    metadata['num_experiments'] = job_context['experiments'].count()
    metadata['quant_sf_only'] = job_context['dataset'].quant_sf_only

    if not job_context['dataset'].quant_sf_only:
        metadata['aggregate_by'] = job_context["dataset"].aggregate_by
        metadata['scale_by'] = job_context["dataset"].scale_by
        # https://github.com/AlexsLemonade/refinebio/pull/421#discussion_r203799646
        # TODO: do something with these.
        # metadata['non_aggregated_files'] = job_context["unsmashable_files"]
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
        sample_accessions = experiment.samples.all().values_list('accession_code', flat=True)
        exp_dict['sample_accession_codes'] = [v for v in sample_accessions]
        experiments[experiment.accession_code] = exp_dict

    metadata['experiments'] = experiments

    return metadata


def write_non_data_files(job_context: Dict) -> Dict:
    """Writes the files that are not the actual data of the dataset.

    This include LICENSE.txt and README.md files and the metadata.

    Adds the key `metadata` to job_context and populates it with all
    the metadata that needs to be written.
    """
    job_context['metadata'] = compile_metadata(job_context)

    shutil.copy("README_DATASET.md", job_context["output_dir"] + "README.md")
    shutil.copy("LICENSE_DATASET.txt", job_context["output_dir"] + "LICENSE.TXT")

    # Write samples metadata to TSV
    try:
        tsv_paths = write_tsv_json(job_context)
        job_context['metadata_tsv_paths'] = tsv_paths
        # Metadata to JSON
        job_context['metadata']['created_at'] = timezone.now().strftime('%Y-%m-%dT%H:%M:%S')
        with open(job_context["output_dir"] + 'aggregated_metadata.json',
                  'w',
                  encoding='utf-8') as metadata_file:
            json.dump(job_context['metadata'], metadata_file, indent=4, sort_keys=True)
    except Exception as e:
        logger.exception("Failed to write metadata TSV!",
                         job_id=job_context['job'].id)
        job_context['metadata_tsv_paths'] = None

    return job_context


def get_experiment_accession(sample_accession_code, dataset_data):
    for experiment_accession, samples in dataset_data.items():
        if sample_accession_code in samples:
            return experiment_accession
    return ""  # Should never happen, because the sample is by definition in the dataset


def _add_annotation_column(annotation_columns, column_name):
    """Add annotation column names in place.
    Any column_name that starts with "refinebio_" will be skipped.
    """

    if not column_name.startswith("refinebio_"):
        annotation_columns.add(column_name)


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


def get_tsv_row_data(sample_metadata, dataset_data):
    """Returns field values based on input sample_metadata.

    Some annotation fields are treated specially because they are more
    important.  See `get_tsv_columns` function above for details.
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
                if sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                   and annotation_key == "characteristic":
                    for pair_dict in annotation_value:
                        if 'category' in pair_dict and 'value' in pair_dict:
                            col_name, col_value = pair_dict['category'], pair_dict['value']
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                # "variable" in ArrayExpress annotation
                elif sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                     and annotation_key == "variable":
                    for pair_dict in annotation_value:
                        if 'name' in pair_dict and 'value' in pair_dict:
                            col_name, col_value = pair_dict['name'], pair_dict['value']
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                # Skip "source" field ArrayExpress sample's annotation
                elif sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                     and annotation_key == "source":
                    continue
                # "characteristics_ch1" in GEO annotation
                elif sample_metadata.get('refinebio_source_database', '') == "GEO" \
                     and annotation_key == "characteristics_ch1":  # array of strings
                    for pair_str in annotation_value:
                        if ':' in pair_str:
                            col_name, col_value = pair_str.split(':', 1)
                            col_value = col_value.strip()
                            _add_annotation_value(row_data, col_name, col_value,
                                                  sample_accession_code)
                # If annotation_value includes only a 'name' key, extract its value directly:
                elif isinstance(annotation_value, dict) \
                     and len(annotation_value) == 1 and 'name' in annotation_value:
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

    row_data["experiment_accession"] = get_experiment_accession(sample_accession_code,
                                                                dataset_data)

    return row_data


def get_tsv_columns(job_context, samples_metadata):
    """Returns an array of strings that will be written as a TSV file's
    header. The columns are based on fields found in samples_metadata.

    Some nested annotation fields are taken out as separate columns
    because they are more important than the others.
    """
    tsv_start = log_state("start get tsv columns", job_context["job"].id)
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
                    if sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                       and annotation_key == "characteristic":
                        for pair_dict in annotation_value:
                            if 'category' in pair_dict and 'value' in pair_dict:
                                _add_annotation_column(annotation_columns, pair_dict['category'])
                    # For ArrayExpress samples, also take out the fields
                    # nested in "variable" as separate columns.
                    elif sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                         and annotation_key == "variable":
                        for pair_dict in annotation_value:
                            if 'name' in pair_dict and 'value' in pair_dict:
                                _add_annotation_column(annotation_columns, pair_dict['name'])
                    # For ArrayExpress samples, skip "source" field
                    elif sample_metadata.get('refinebio_source_database', '') == "ARRAY_EXPRESS" \
                         and annotation_key == "source":
                        continue
                    # For GEO samples, take out the fields nested in
                    # "characteristics_ch1" as separate columns.
                    elif sample_metadata.get('refinebio_source_database', '') == "GEO" \
                         and annotation_key == "characteristics_ch1":  # array of strings
                        for pair_str in annotation_value:
                            if ':' in pair_str:
                                tokens = pair_str.split(':', 1)
                                _add_annotation_column(annotation_columns, tokens[0])
                    # Saves all other annotation fields in separate columns
                    else:
                        _add_annotation_column(annotation_columns, annotation_key)

    # Return sorted columns, in which "refinebio_accession_code" and "experiment_accession" are
    # always first, followed by the other refinebio columns (in alphabetic order), and
    # annotation columns (in alphabetic order) at the end.
    refinebio_columns.discard('refinebio_accession_code')
    log_state("end get tsv columns", job_context["job"].id, tsv_start)
    return ['refinebio_accession_code', 'experiment_accession'] + sorted(refinebio_columns) \
        + sorted(annotation_columns)


def write_tsv_json(job_context):
    """Writes tsv files on disk.
    If the dataset is aggregated by species, also write species-level
    JSON file.
    """
    # Avoid pulling this out of job_context repeatedly.
    metadata = job_context['metadata']

    # Uniform TSV header per dataset
    columns = get_tsv_columns(job_context, metadata['samples'])

    # Per-Experiment Metadata
    if job_context["dataset"].aggregate_by == "EXPERIMENT":
        tsv_paths = []
        for experiment_title, experiment_data in metadata['experiments'].items():
            experiment_dir = job_context["output_dir"] + experiment_title + '/'
            experiment_dir = experiment_dir.encode('ascii', 'ignore')
            os.makedirs(experiment_dir, exist_ok=True)
            tsv_path = experiment_dir.decode("utf-8") + 'metadata_' + experiment_title + '.tsv'
            tsv_path = tsv_path.encode('ascii', 'ignore')
            tsv_paths.append(tsv_path)
            with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t', extrasaction='ignore')
                dw.writeheader()
                for sample_accession_code, sample_metadata in metadata['samples'].items():
                    if sample_accession_code in experiment_data['sample_accession_codes']:
                        row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
                        dw.writerow(row_data)
        return tsv_paths
    # Per-Species Metadata
    elif job_context["dataset"].aggregate_by == "SPECIES":
        tsv_paths = []
        for species in job_context['group_by_keys']:
            species_dir = job_context["output_dir"] + species + '/'
            os.makedirs(species_dir, exist_ok=True)
            samples_in_species = []
            tsv_path = species_dir + "metadata_" + species + '.tsv'
            tsv_paths.append(tsv_path)
            with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
                # See http://www.lucainvernizzi.net/blog/2015/08/03/8x-speed-up-for-python-s-csv-dictwriter/
                # about extrasaction.
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t', extrasaction='ignore')
                dw.writeheader()
                i = 0
                for sample_metadata in metadata['samples'].values():
                    if sample_metadata.get('refinebio_organism', '') == species:
                        row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
                        dw.writerow(row_data)
                        samples_in_species.append(sample_metadata)

                    i = i + 1
                    if i % 1000 == 0:
                        progress_template = ('Done with {0} out of {1} lines of metadata '
                                             'for species {2}')
                        log_state(progress_template.format(i, len(metadata['samples']), species),
                                  job_context['job'].id)

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
        all_dir = job_context["output_dir"] + "ALL/"
        os.makedirs(all_dir, exist_ok=True)
        tsv_path = all_dir + 'metadata_ALL.tsv'
        with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
            dw = csv.DictWriter(tsv_file, columns, delimiter='\t', extrasaction='ignore')
            dw.writeheader()
            for sample_metadata in metadata['samples'].values():
                row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
                dw.writerow(row_data)
        return [tsv_path]


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


def sync_quant_files(output_path, samples: List[Sample]):
    """ Takes a list of ComputedFiles and copies the ones that are quant files to the provided directory.
        Returns the total number of samples that were included """
    num_samples = 0

    page_size = 100
    # split the samples in groups and download each one individually
    pool = multiprocessing.Pool(processes=MULTIPROCESSING_WORKER_COUNT)

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

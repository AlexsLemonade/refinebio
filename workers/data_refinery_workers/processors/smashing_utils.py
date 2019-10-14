# -*- coding: utf-8 -*-

import csv
import itertools
import logging
import multiprocessing
import os
import psutil
import rpy2.robjects as ro
import shutil
import simplejson as json
import time

from django.utils import timezone
from pathlib import Path
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as rlang
from rpy2.robjects.packages import importr
from typing import Dict, List
import numpy as np
import pandas as pd

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputedFile
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils


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

MULTIPROCESSING_WORKER_COUNT = max(1, multiprocessing.cpu_count()/2 - 1)
MULTIPROCESSING_CHUNK_SIZE = 2000


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
    #       fgenesh2_kg.7__3016__AT5G35080.1 (via http://plants.ensembl.org/Arabidopsis_lyrata/ \
    #       Gene/Summary?g=fgenesh2_kg.7__3016__AT5G35080.1;r=7:17949732-17952000;t=fgenesh2_kg. \
    #       7__3016__AT5G35080.1;db=core)
    data.index = data.index.str.replace(r"(\.[^.]*)$", '')

    # Squish duplicated rows together.
    # XXX/TODO: Is mean the appropriate method here?
    #           We can make this an option in future.
    # Discussion here: https://github.com/AlexsLemonade/refinebio/issues/186#issuecomment-395516419
    data = data.groupby(data.index, sort=False).mean()

    return data


def process_frame(inputs) -> Dict:
    (
        work_dir,
        computed_file,
        sample_accession_code,
        dataset_id,
        aggregate_by,
        index,
        job_id
    ) = inputs

    logger.debug('processing frame {}'.format(index), job_id=job_id)

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

        data = _load_and_sanitize_file(computed_file_path)

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
    finally:
        # Delete before archiving the work dir
        if computed_file_path and os.path.exists(computed_file_path):
            os.remove(computed_file_path)

    return smashable(data)


def process_frames_for_key(key: str, input_files: List[ComputedFile], job_context: Dict) -> Dict:
    job_context['original_merged'] = pd.DataFrame()

    start_frames = log_state("building frames for species or experiment {}".format(key),
                             job_context["job"])

    def get_frame_inputs():
        """Helper method to create a generator object that when iterated on returns
        a tuple for smashing_utils.process_frame method to process a single frame.
        It returns a tuple becuase process_frame is called with multiprocessing.Pool.map
        which restricts how arguments are passed"""
        for index, (computed_file, sample) in enumerate(input_files):
            # Don't pass job_context to worker threads because
            # for some reason it causes them to open database
            # connections. Not yet sure why...
            yield (
                job_context["work_dir"],
                computed_file,
                sample.accession_code,
                job_context['dataset'].id,
                job_context['dataset'].aggregate_by,
                index,
                job_context["job"].id
            )

    def get_chunked_frame_inputs():
        """Helper method for chunking the generator results returned from get_frame_inputs.
        We want to be able to pass a single generator that returns a chunk of
        frame_inputs to a multiprocessing.Pool.map and we want a generator that will
        return those generators. This function returns a generator that returns the chunked
        generator when iterated on which returns a single frame_input tuple when iterated on.
        Essentially it returns a lazily loaded list of lists of tuples.
        """
        # get the resulting generator that returns all tuples for
        # smashing_utils.process_frame
        source = get_frame_inputs()
        # get a generator for the first chunk - at this point we might only have one
        iterator = itertools.islice(source, MULTIPROCESSING_CHUNK_SIZE)
        while True:
            # we need to look and see if there is a value in the chunk
            # if there is nothing return None and we can exit this loop
            first_item_in_next_chunk = next(iterator, None)
            if first_item_in_next_chunk is None:
                return
            # at this point we know that our chunk has at least one item in it
            # we want to add back the result that we took out and reset what we had
            source, iterator = itertools.tee(itertools.chain(
                [first_item_in_next_chunk],
                source
            ))
            # yield a generator for the current chunk which terminates
            # at MULTIPROCESSING_CHUNK_SIZE
            yield itertools.islice(iterator, MULTIPROCESSING_CHUNK_SIZE)
            # we want to move the source generator to start where
            # the last chunk left off
            next(itertools.islice(
                 source,
                 MULTIPROCESSING_CHUNK_SIZE,
                 MULTIPROCESSING_CHUNK_SIZE+1
            ), None)

    # here we create a generator that as it consumes a chunk generator
    # it passes that to a worker pool and yields the processed frames
    def get_processed_frames_with_pool(pool):
        for chunk_of_frames in get_chunked_frame_inputs():
            yield pool.map(process_frame, chunk_of_frames)

    # Build up a list of microarray frames and a list of
    # rnaseq frames and then combine them so they're sorted
    # out.
    job_context['microarray_frames'] = []
    job_context['rnaseq_frames'] = []

    # take one fewer than 1/2 the total available threads
    # also make the minimum threads 1
    worker_pool = multiprocessing.Pool(processes=MULTIPROCESSING_WORKER_COUNT)

    chunks_from_pool = get_processed_frames_with_pool(worker_pool)

    # use the results from each frame from each chunk
    # we are using itertool.chain.from_iterable as a
    # way to flatten the results
    for frame in itertools.chain.from_iterable(chunks_from_pool):
        if frame['technology'] == 'microarray':
            job_context['microarray_frames'].append(frame['dataframe'])
        elif frame['technology'] == 'rnaseq':
            job_context['rnaseq_frames'].append(frame['dataframe'])

        if frame['unsmashable']:
            job_context['unsmashable_files'].append(frame['unsmashable_file'])

    # clean up the pool when we are done
    worker_pool.close()
    worker_pool.join()

    job_context['num_samples'] = job_context['num_samples'] \
                                 + len(job_context['microarray_frames'])
    job_context['num_samples'] = job_context['num_samples'] + len(job_context['rnaseq_frames'])

    log_state("set frames for key {}".format(key), job_context["job"], start_frames)

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

    Expects the key `metadata` in job_context to be populated with all
    the metadata that needs to be written.
    """
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
    tsv_start = log_state("start get tsv columns", job_context["job"])
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
    log_state("end get tsv columns", job_context["job"], tsv_start)
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
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
                dw.writeheader()
                for sample_accession_code, sample_metadata in metadata['samples'].items():
                    if sample_accession_code in experiment_data['sample_accession_codes']:
                        row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
                        dw.writerow(row_data)
        return tsv_paths
    # Per-Species Metadata
    elif job_context["dataset"].aggregate_by == "SPECIES":
        tsv_paths = []
        for species in job_context['input_files'].keys():
            species_dir = job_context["output_dir"] + species + '/'
            os.makedirs(species_dir, exist_ok=True)
            samples_in_species = []
            tsv_path = species_dir + "metadata_" + species + '.tsv'
            tsv_paths.append(tsv_path)
            with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
                dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
                dw.writeheader()
                for sample_metadata in metadata['samples'].values():
                    if sample_metadata.get('refinebio_organism', '') == species:
                        row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
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
        all_dir = job_context["output_dir"] + "ALL/"
        os.makedirs(all_dir, exist_ok=True)
        tsv_path = all_dir + 'metadata_ALL.tsv'
        with open(tsv_path, 'w', encoding='utf-8') as tsv_file:
            dw = csv.DictWriter(tsv_file, columns, delimiter='\t')
            dw.writeheader()
            for sample_metadata in metadata['samples'].values():
                row_data = get_tsv_row_data(sample_metadata, job_context["dataset"].data)
                dw.writerow(row_data)
        return [tsv_path]

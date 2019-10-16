import logging
import os
import psutil
import shutil
import time

import numpy as np
import pandas as pd
from fancyimpute import IterativeSVD
import simplejson as json

from datetime import datetime
from django.utils import timezone
from typing import Dict

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Pipeline,
    Organism
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils, smashing_utils


pd.set_option('mode.chained_assignment', None)


S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
BYTES_IN_GB = 1024 * 1024 * 1024
SMASHING_DIR = "/home/user/data_store/smashed/"
logger = get_and_configure_logger(__name__)
### DEBUG ###
logger.setLevel(logging.getLevelName('DEBUG'))


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


def _prepare_input(job_context: Dict) -> Dict:
    start_time = log_state("prepare input", job_context["job"])

    job_context = smashing_utils.prepare_files(job_context)
    if job_context['job'].success is False:
        smashing_utils.log_failure(job_context, "Unable to run smashing_utils.prepare_files.")
        return job_context

    # Compendia jobs only run for one organism, so we know the only
    # key will be the organism name, unless of course we've already failed.
    if job_context['job'].success is not False:
        job_context["organism_name"] = list(job_context['input_files'].keys())[0]

    log_state("prepare input done", job_context["job"], start_time)
    return job_context


def _prepare_frames(job_context: Dict) -> Dict:
    start_prepare_frames = log_state("start _prepare_frames", job_context["job"])

    try:
        job_context['unsmashable_files'] = []
        job_context['num_samples'] = 0

        # Smash all of the sample sets
        logger.debug("About to smash!",
                     dataset_count=len(job_context['dataset'].data),
                     job_id=job_context['job'].id)

        # Once again, `key` is either a species name or an experiment accession
        for key, input_files in job_context['input_files'].items():
            job_context = smashing_utils.process_frames_for_key(key, input_files, job_context)
            # if len(job_context['all_frames']) < 1:
            # TODO: Enable this check?

    except Exception as e:
        logger.exception("Could not prepare frames for compendia.",
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

    log_state("end _prepare_frames", job_context["job"], start_prepare_frames)
    return job_context


def _perform_imputation(job_context: Dict) -> Dict:
    """

    Take the inputs and perform the primary imputation.

    Via https://github.com/AlexsLemonade/refinebio/issues/508#issuecomment-435879283:
     - Combine all microarray samples with a full join to form a
       microarray_expression_matrix (this may end up being a DataFrame).
     - Combine all RNA-seq samples (lengthScaledTPM) with a full outer join
       to form a rnaseq_expression_matrix.
     - Calculate the sum of the lengthScaledTPM values for each row (gene) of
       the rnaseq_expression_matrix (rnaseq_row_sums).
     - Calculate the 10th percentile of rnaseq_row_sums
     - Drop all rows in rnaseq_expression_matrix with a row sum < 10th percentile of
       rnaseq_row_sums; this is now filtered_rnaseq_matrix
     - log2(x + 1) transform filtered_rnaseq_matrix; this is now log2_rnaseq_matrix
     - Set all zero values in log2_rnaseq_matrix to NA, but make sure to keep track of
       where these zeroes are
     - Perform a full outer join of microarray_expression_matrix and
       log2_rnaseq_matrix; combined_matrix
     - Remove genes (rows) with >30% missing values in combined_matrix
     - Remove samples (columns) with >50% missing values in combined_matrix
     - "Reset" zero values that were set to NA in RNA-seq samples (i.e., make these zero
       again) in combined_matrix
     - Transpose combined_matrix; transposed_matrix
     - Perform imputation of missing values with IterativeSVD (rank=10) on
       the transposed_matrix; imputed_matrix
        -- with specified svd algorithm or skip
     - Untranspose imputed_matrix (genes are now rows, samples are now columns)
     - Quantile normalize imputed_matrix where genes are rows and samples are columns

    """
    imputation_start = log_state("start perform imputation", job_context["job"])
    job_context['time_start'] = timezone.now()

    microarray_start = log_state("start microarray concatenation", job_context["job"])

    # Combine all microarray samples with a full outer join to form a
    # microarray_expression_matrix (a DataFrame).
    microarray_expression_matrix = pd.concat(job_context.pop('microarray_frames'),
                                             axis=1,
                                             keys=None,
                                             join='outer',
                                             copy=False,
                                             sort=True)

    log_state("end microarray concatenation", job_context["job"], microarray_start)
    rnaseq_start = log_state("start rnaseq concatenation", job_context["job"])

    # Combine all RNA-seq samples (lengthScaledTPM) with a full outer
    # join to form a rnaseq_expression_matrix (a DataFrame).
    rnaseq_expression_matrix = pd.concat(job_context.pop('rnaseq_frames'),
                                         axis=1,
                                         keys=None,
                                         join='outer',
                                         copy=False,
                                         sort=True)

    log_state("end rnaseq concatenation", job_context["job"], rnaseq_start)
    rnaseq_row_sums_start = log_state("start rnaseq row sums", job_context["job"])

    # Calculate the sum of the lengthScaledTPM values for each row
    # (gene) of the rnaseq_expression_matrix (rnaseq_row_sums)
    rnaseq_row_sums = np.sum(rnaseq_expression_matrix, axis=1)

    log_state("end rnaseq row sums", job_context["job"], rnaseq_row_sums_start)
    rnaseq_decile_start = log_state("start rnaseq decile", job_context["job"])

    # Calculate the 10th percentile of rnaseq_row_sums
    rnaseq_tenth_percentile = np.percentile(rnaseq_row_sums, 10)

    log_state("end rnaseq decile", job_context["job"], rnaseq_decile_start)
    drop_start = log_state("drop all rows", job_context["job"])
    # Drop all rows in rnaseq_expression_matrix with a row sum < 10th
    # percentile of rnaseq_row_sums; this is now
    # filtered_rnaseq_matrix
    # TODO: This is probably a better way to do this with `np.where`
    rows_to_filter = []
    for (x, sum_val) in rnaseq_row_sums.items():
        if sum_val < rnaseq_tenth_percentile:
            rows_to_filter.append(x)

    del rnaseq_row_sums

    log_state("actually calling drop()", job_context["job"])

    filtered_rnaseq_matrix = rnaseq_expression_matrix.drop(rows_to_filter)

    log_state("end drop all rows", job_context["job"], drop_start)
    log2_start = log_state("start log2", job_context["job"])

    # log2(x + 1) transform filtered_rnaseq_matrix; this is now log2_rnaseq_matrix
    filtered_rnaseq_matrix_plus_one = filtered_rnaseq_matrix + 1
    log2_rnaseq_matrix = np.log2(filtered_rnaseq_matrix_plus_one)
    del filtered_rnaseq_matrix_plus_one
    del filtered_rnaseq_matrix

    log_state("end log2", job_context["job"], log2_start)
    cache_start = log_state("start caching zeroes", job_context["job"])

    # Cache our RNA-Seq zero values
    cached_zeroes = {}
    for column in log2_rnaseq_matrix.columns:
        cached_zeroes[column] = log2_rnaseq_matrix.index[np.where(log2_rnaseq_matrix[column] == 0)]

    log_state("end caching zeroes", job_context["job"], cache_start)
    outer_merge_start = log_state("start outer merge", job_context["job"])

    # Set all zero values in log2_rnaseq_matrix to NA, but make sure
    # to keep track of where these zeroes are
    log2_rnaseq_matrix[log2_rnaseq_matrix == 0] = np.nan

    # Perform a full outer join of microarray_expression_matrix and
    # log2_rnaseq_matrix; combined_matrix
    combined_matrix = microarray_expression_matrix.merge(log2_rnaseq_matrix,
                                                         how='outer',
                                                         left_index=True,
                                                         right_index=True)
    del microarray_expression_matrix

    log_state("end outer merge", job_context["job"], outer_merge_start)
    drop_na_genes_start = log_state("start drop NA genes", job_context["job"])

    # # Visualize Prefiltered
    # output_path = job_context['output_dir'] + "pre_filtered_" + str(time.time()) + ".png"
    # visualized_prefilter = visualize.visualize(combined_matrix.copy(), output_path)

    # Remove genes (rows) with <=70% present values in combined_matrix
    thresh = combined_matrix.shape[1] * .7  # (Rows, Columns)
    # Everything below `thresh` is dropped
    row_filtered_matrix = combined_matrix.dropna(axis='index',
                                                 thresh=thresh)
    del combined_matrix
    del thresh

    log_state("end drop NA genes", job_context["job"], drop_na_genes_start)
    drop_na_samples_start = log_state("start drop NA samples", job_context["job"])

    # # Visualize Row Filtered
    # output_path = job_context['output_dir'] + "row_filtered_" + str(time.time()) + ".png"
    # visualized_rowfilter = visualize.visualize(row_filtered_matrix.copy(), output_path)

    # Remove samples (columns) with <50% present values in combined_matrix
    # XXX: Find better test data for this!
    col_thresh = row_filtered_matrix.shape[0] * .5
    row_col_filtered_matrix_samples = row_filtered_matrix.dropna(axis='columns',
                                                                 thresh=col_thresh)
    row_col_filtered_matrix_samples_index = row_col_filtered_matrix_samples.index
    row_col_filtered_matrix_samples_columns = row_col_filtered_matrix_samples.columns

    log_state("end drop NA genes", job_context["job"], drop_na_samples_start)
    replace_zeroes_start = log_state("start replace zeroes", job_context["job"])

    del row_filtered_matrix

    # # Visualize Row and Column Filtered
    # output_path = job_context['output_dir'] + "row_col_filtered_" + str(time.time()) + ".png"
    # visualized_rowcolfilter = visualize.visualize(row_col_filtered_matrix_samples.copy(),
    #                                               output_path)

    # "Reset" zero values that were set to NA in RNA-seq samples
    # (i.e., make these zero again) in combined_matrix
    for column in cached_zeroes.keys():
        zeroes = cached_zeroes[column]

        # Skip purged columns
        if column not in row_col_filtered_matrix_samples:
            continue

        # Place the zero
        try:
            # This generates a warning, so use loc[] instead
            # row_col_filtered_matrix_samples[column].replace(zeroes, 0.0, inplace=True)
            zeroes_list = zeroes.tolist()
            new_index_list = row_col_filtered_matrix_samples_index.tolist()
            new_zeroes = list(set(new_index_list) & set(zeroes_list))
            row_col_filtered_matrix_samples[column].loc[new_zeroes] = 0.0
        except Exception as e:
            logger.warn("Error when replacing zero")
            continue

    log_state("end replace zeroes", job_context["job"], replace_zeroes_start)
    transposed_zeroes_start = log_state("start replacing transposed zeroes", job_context["job"])

    # Label our new replaced data
    combined_matrix_zero = row_col_filtered_matrix_samples
    del row_col_filtered_matrix_samples

    transposed_matrix_with_zeros = combined_matrix_zero.T
    del combined_matrix_zero

    # Remove -inf and inf
    # This should never happen, but make sure it doesn't!
    transposed_matrix = transposed_matrix_with_zeros.replace([np.inf, -np.inf], np.nan)
    del transposed_matrix_with_zeros

    log_state("end replacing transposed zeroes", job_context["job"], transposed_zeroes_start)

    # Store the absolute/percentages of imputed values
    matrix_sum = transposed_matrix.isnull().sum()
    percent = (matrix_sum / transposed_matrix.isnull().count()).sort_values(ascending=False)
    total_percent_imputed = sum(percent) / len(transposed_matrix.count())
    job_context['total_percent_imputed'] = total_percent_imputed
    logger.info("Total percentage of data to impute!", total_percent_imputed=total_percent_imputed)

    # Perform imputation of missing values with IterativeSVD (rank=10) on the
    # transposed_matrix; imputed_matrix
    svd_algorithm = job_context['dataset'].svd_algorithm
    if svd_algorithm != 'NONE':
        svd_start = log_state("start SVD", job_context["job"])

        logger.info("IterativeSVD algorithm: %s" % svd_algorithm)
        svd_algorithm = str.lower(svd_algorithm)
        imputed_matrix = IterativeSVD(
            rank=10,
            svd_algorithm=svd_algorithm
        ).fit_transform(transposed_matrix)

        svd_start = log_state("end SVD", job_context["job"], svd_start)
    else:
        imputed_matrix = transposed_matrix
        logger.info("Skipping IterativeSVD")
    del transposed_matrix

    untranspose_start = log_state("start untranspose", job_context["job"])

    # Untranspose imputed_matrix (genes are now rows, samples are now columns)
    untransposed_imputed_matrix = imputed_matrix.T
    del imputed_matrix

    # Convert back to Pandas
    untransposed_imputed_matrix_df = pd.DataFrame.from_records(untransposed_imputed_matrix)
    untransposed_imputed_matrix_df.index = row_col_filtered_matrix_samples_index
    untransposed_imputed_matrix_df.columns = row_col_filtered_matrix_samples_columns
    del untransposed_imputed_matrix
    del row_col_filtered_matrix_samples_index
    del row_col_filtered_matrix_samples_columns
    # Quantile normalize imputed_matrix where genes are rows and samples are columns
    job_context['organism'] = Organism.get_object_for_name(job_context['organism_name'])
    job_context['merged_no_qn'] = untransposed_imputed_matrix_df
    # output_path = job_context['output_dir'] + "compendia_no_qn_" + str(time.time()) + ".png"
    # visualized_merged_no_qn = visualize.visualize(untransposed_imputed_matrix_df.copy(),
    #                                               output_path)

    log_state("end untranspose", job_context["job"], untranspose_start)
    quantile_start = log_state("start quantile normalize", job_context["job"])

    # Perform the Quantile Normalization
    job_context = smashing_utils.quantile_normalize(job_context, ks_check=False)

    log_state("end quantile normalize", job_context["job"], quantile_start)

    # Visualize Final Compendia
    # output_path = job_context['output_dir'] + "compendia_with_qn_" + str(time.time()) + ".png"
    # visualized_merged_qn = visualize.visualize(job_context['merged_qn'].copy(), output_path)

    job_context['time_end'] = timezone.now()
    job_context['formatted_command'] = "create_compendia.py"
    log_state("end prepare imputation", job_context["job"], imputation_start)
    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    """
    Store and host the result as a ComputationalResult object.
    """
    result_start = log_state("start create result object", job_context["job"])
    result = ComputationalResult()
    result.commands.append(" ".join(job_context['formatted_command']))
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context['time_start']
    result.time_end = job_context['time_end']
    try:
        processor_key = "CREATE_COMPENDIA"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    # Write the compendia dataframe to a file
    job_context['csv_outfile'] = job_context['output_dir'] + job_context['organism_name'] + '.tsv'
    job_context['merged_qn'].to_csv(job_context['csv_outfile'], sep='\t', encoding='utf-8')
    compendia_tsv_computed_file = ComputedFile()
    compendia_tsv_computed_file.absolute_file_path = job_context['csv_outfile']
    compendia_tsv_computed_file.filename = job_context['csv_outfile'].split('/')[-1]
    compendia_tsv_computed_file.calculate_sha1()
    compendia_tsv_computed_file.calculate_size()
    compendia_tsv_computed_file.is_smashable = False
    compendia_tsv_computed_file.is_qn_target = False
    compendia_tsv_computed_file.result = result
    compendia_tsv_computed_file.save()

    organism_key = list(job_context['samples'].keys())[0]
    annotation = ComputationalResultAnnotation()
    annotation.result = result

    annotation.data = {
        "organism_id": job_context['samples'][organism_key][0].organism_id,
        "organism_name": job_context['organism_name'],
        "is_qn": False,
        "is_compendia": True,
        "samples": [sample.accession_code for sample in job_context["samples"][organism_key]],
        "num_samples": len(job_context["samples"][organism_key]),
        "experiment_accessions": [e.accession_code for e in job_context['experiments']],
        "total_percent_imputed": job_context['total_percent_imputed']
    }
    annotation.save()

    # Save the related metadata file
    metadata_computed_file = ComputedFile()
    metadata_computed_file.absolute_file_path = job_context['metadata_tsv_paths'][0]
    metadata_computed_file.filename = job_context['metadata_tsv_paths'][0].split('/')[-1]
    metadata_computed_file.calculate_sha1()
    metadata_computed_file.calculate_size()
    metadata_computed_file.is_smashable = False
    metadata_computed_file.is_qn_target = False
    metadata_computed_file.result = result
    metadata_computed_file.save()

    # Create the resulting archive
    final_zip_base = SMASHING_DIR + str(job_context["dataset"].pk) + "_compendia"
    # Copy LICENSE.txt and correct README.md files.
    if job_context["dataset"].quant_sf_only:
        readme_file = "/home/user/README_QUANT.md"
    else:
        readme_file = "/home/user/README_NORMALIZED.md"

    shutil.copy(readme_file, job_context["output_dir"] + "/README.md")
    shutil.copy("/home/user/LICENSE_DATASET.txt", job_context["output_dir"] + "/LICENSE.TXT")
    archive_path = shutil.make_archive(final_zip_base, 'zip', job_context["output_dir"])

    # Save the related metadata file
    organism = job_context['samples'][organism_key][0].organism

    try:
        last_compendia = ComputedFile.objects.filter(
                                    is_compendia=True,
                                    compendia_organism=organism).order_by('-compendia_version')[-1]
        compendia_version = last_compendia.compendia_version + 1
    except Exception as e:
        # This is the first compendia for this Organism
        compendia_version = 1

    archive_computed_file = ComputedFile()
    archive_computed_file.absolute_file_path = archive_path
    archive_computed_file.filename = archive_path.split('/')[-1]
    archive_computed_file.calculate_sha1()
    archive_computed_file.calculate_size()
    archive_computed_file.is_smashable = False
    archive_computed_file.is_qn_target = False
    archive_computed_file.result = result
    archive_computed_file.is_compendia = True
    archive_computed_file.quant_sf_only = job_context["dataset"].quant_sf_only
    archive_computed_file.svd_algorithm = job_context["dataset"].svd_algorithm
    archive_computed_file.compendia_organism = job_context['samples'][organism_key][0].organism
    archive_computed_file.compendia_version = compendia_version
    archive_computed_file.save()

    logger.info("Compendia created!",
                archive_path=archive_path,
                organism_name=job_context['organism_name'])

    # Upload the result to S3
    timestamp = str(int(time.time()))
    key = job_context['organism_name'] + "_" + str(compendia_version) + "_" + timestamp + ".zip"
    archive_computed_file.sync_to_s3(S3_BUCKET_NAME, key)

    job_context['result'] = result
    job_context['computed_files'] = [compendia_tsv_computed_file,
                                     metadata_computed_file,
                                     archive_computed_file]
    job_context['success'] = True

    log_state("end create result object", job_context["job"], result_start)

    return job_context


def create_compendia(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.CREATE_COMPENDIA.value)
    job_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                                     [utils.start_job,
                                      _prepare_input,
                                      _prepare_frames,
                                      _perform_imputation,
                                      smashing_utils.write_non_data_files,
                                      _create_result_objects,
                                      utils.end_job])
    return job_context

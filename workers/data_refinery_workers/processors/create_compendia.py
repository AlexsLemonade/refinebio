import os
import random
import string
import subprocess
import time
import warnings

import numpy as np
import pandas as pd
from fancyimpute import KNN, BiScaler, SoftImpute, IterativeSVD

from django.utils import timezone
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Pipeline,
    Processor,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    Organism
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils, smasher


S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
logger = get_and_configure_logger(__name__)


def _prepare_input(job_context: Dict) -> Dict:

    # We're going to use the smasher outside of the smasher.
    # I'm not crazy about this yet. Maybe refactor later,
    # but I need the data now.
    job_context = smasher._prepare_files(job_context)
    job_context = smasher._smash(job_context)

    if not 'final_frame' in job_context.keys():
        logger.error("Unable to prepare files for creating compendia.",
            job_id=job_context['job'].id)
        job_context["job"].failure_reason = "Couldn't prepare files creating compendia."
        job_context['success'] = False
        return job_context

    # Prep the two data types for imputation
    og_merged = job_context['original_merged']
    job_context['microarray_inputs'] = og_merged.copy()
    job_context['rnaseq_inputs'] = og_merged.copy()

    for column_name in og_merged.columns:
        if column_name in job_context['technologies']['microarray']:
            del job_context['rnaseq_inputs'][column_name]
        else:
            del job_context['microarray_inputs'][column_name]

    # work_dir is already created by smasher._prepare_files
    outfile_base = job_context['work_dir'] + str(time.time()).split('.')[0]
    outfile = outfile_base + '.tsv'
    job_context['final_frame'].to_csv(outfile, sep='\t', encoding='utf-8')
    job_context['smashed_file'] = outfile
    job_context['target_file'] = outfile_base + '_target.tsv'

    return job_context

def _perform_imputation(job_context: Dict) -> Dict:
    """

    Take the inputs and perform the primary imputation.

    Via https://github.com/AlexsLemonade/refinebio/issues/508#issuecomment-435879283: 
     - Combine all microarray samples with a full join to form a microarray_expression_matrix (this may end up being a DataFrame)
     - Combine all RNA-seq samples (lengthScaledTPM) with a full outer join to form a rnaseq_expression_matrix
     - Calculate the sum of the lengthScaledTPM values for each row (gene) of the rnaseq_expression_matrix (rnaseq_row_sums)
     - Calculate the 10th percentile of rnaseq_row_sums
     - Drop all rows in rnaseq_expression_matrix with a row sum < 10th percentile of rnaseq_row_sums; this is now filtered_rnaseq_matrix
     - log2(x + 1) transform filtered_rnaseq_matrix; this is now log2_rnaseq_matrix
     - Set all zero values in log2_rnaseq_matrix to NA, but make sure to keep track of where these zeroes are
     - Perform a full outer join of microarray_expression_matrix and log2_rnaseq_matrix; combined_matrix
     - Remove genes (rows) with >30% missing values in combined_matrix
     - Remove samples (columns) with >50% missing values in combined_matrix
     - "Reset" zero values that were set to NA in RNA-seq samples (i.e., make these zero again) in combined_matrix
     - Transpose combined_matrix; transposed_matrix
     - Perform imputation of missing values with IterativeSVD (rank=10) on the transposed_matrix; imputed_matrix
     - Untranspose imputed_matrix (genes are now rows, samples are now columns)
     - Quantile normalize imputed_matrix where genes are rows and samples are columns

    """
    job_context['time_start'] = timezone.now()

    # Combine all microarray samples with a full join to form a microarray_expression_matrix (this may end up being a DataFrame)
    microarray_expression_matrix = job_context['microarray_inputs']

    # Combine all RNA-seq samples (lengthScaledTPM) with a full outer join to form a rnaseq_expression_matrix
    rnaseq_expression_matrix = job_context['rnaseq_inputs']

    # Calculate the sum of the lengthScaledTPM values for each row (gene) of the rnaseq_expression_matrix (rnaseq_row_sums)
    rnaseq_row_sums = np.sum(rnaseq_expression_matrix, axis=1)

    # Calculate the 10th percentile of rnaseq_row_sums
    rnaseq_tenth_percentile = np.percentile(rnaseq_row_sums, 10)

    # Drop all rows in rnaseq_expression_matrix with a row sum < 10th percentile of rnaseq_row_sums; this is now filtered_rnaseq_matrix
    # TODO: This is probably a better way to do this with `np.where`
    rows_to_filter = []
    for (x, sum_val) in rnaseq_row_sums.items():
        if sum_val < rnaseq_tenth_percentile:
            rows_to_filter.append(x)

    filtered_rnaseq_matrix = rnaseq_expression_matrix.drop(rows_to_filter)

    # log2(x + 1) transform filtered_rnaseq_matrix; this is now log2_rnaseq_matrix
    filtered_rnaseq_matrix_plus_one = filtered_rnaseq_matrix + 1
    log2_rnaseq_matrix = np.log2(filtered_rnaseq_matrix)

    # Set all zero values in log2_rnaseq_matrix to NA, but make sure to keep track of where these zeroes are
    log2_rnaseq_matrix[log2_rnaseq_matrix==0]=np.nan

    # Perform a full outer join of microarray_expression_matrix and log2_rnaseq_matrix; combined_matrix
    combined_matrix = pd.merge(microarray_expression_matrix, log2_rnaseq_matrix, how='outer', left_index=True, right_index=True)

    # Remove genes (rows) with >30% missing values in combined_matrix
    thresh = len(combined_matrix.columns) * .6
    filtered_combined_matrix = combined_matrix.dropna(axis=0, thresh=thresh, how='any')

    # Remove samples (columns) with >50% missing values in combined_matrix
    # XXX: Find better test data for this!
    thresh = len(combined_matrix.columns) * .5
    filtered_combined_matrix_samples = filtered_combined_matrix.dropna(axis=1, thresh=thresh, how='any')

    # "Reset" zero values that were set to NA in RNA-seq samples (i.e., make these zero again) in combined_matrix
    #combined_matrix_zero = filtered_combined_matrix_samples.fillna(value=0)

    # Transpose combined_matrix; transposed_matrix
    transposed_matrix = filtered_combined_matrix_samples.transpose()

    # Remove -inf and inf
    transposed_matrix = transposed_matrix.replace([np.inf, -np.inf], np.nan)

    # Perform imputation of missing values with IterativeSVD (rank=10) on the transposed_matrix; imputed_matrix
    imputed_matrix = IterativeSVD(rank=10).fit_transform(transposed_matrix)

    # Untranspose imputed_matrix (genes are now rows, samples are now columns)
    untransposed_imputed_matrix = imputed_matrix.transpose()

    # Convert back to Pandas
    untransposed_imputed_matrix_df = pd.DataFrame.from_records(untransposed_imputed_matrix)
    untransposed_imputed_matrix_df.index = filtered_combined_matrix_samples.index
    untransposed_imputed_matrix_df.columns = filtered_combined_matrix_samples.columns

    # Quantile normalize imputed_matrix where genes are rows and samples are columns
    # XXX: Refactor QN target acquisition and application before doing this
    job_context['organism'] = Organism.get_object_for_name(list(job_context['input_files'].keys())[0])
    job_context['merged_no_qn'] = untransposed_imputed_matrix_df
    job_context = smasher._quantile_normalize(job_context)

    job_context['formatted_command'] = "create_compendia.py" # ???

    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    """
    Store and host the result as a ComputationalResult object.
    """

    result = ComputationalResult()
    result.commands.append(" ".join(job_context['formatted_command']))
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context['time_start']
    result.time_end = job_context['time_end']
    try:
        processor_key = "COMPENDIA"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context['target_file']
    computed_file.filename = job_context['target_file'].split('/')[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.is_smashable = False
    computed_file.is_qn_target = True
    computed_file.result = result
    computed_file.save()

    annotation = ComputationalResultAnnotation()
    annotation.result = result
    annotation.data = {
        "organism_id": job_context['samples']['ALL'][0].organism_id,
        "is_qn": False,
        "platform_accession_code": job_context['samples']['ALL'][0].platform_accession_code,
        "samples": [sample.accession_code for sample in job_context["samples"]["ALL"]]
    }
    annotation.save()

    # TODO: upload this to a public read bucket.
    # https://github.com/AlexsLemonade/refinebio/issues/586
    job_context['result'] = result
    job_context['computed_files'] = [computed_file]
    job_context['success'] = True
    return job_context

def create_compendia(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.COMPENDIA.value)
    job_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                       [utils.start_job,
                        _prepare_input,
                        _perform_imputation,
                        _create_result_objects,
                        utils.end_job])
    return job_context

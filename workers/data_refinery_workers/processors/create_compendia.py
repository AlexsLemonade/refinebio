import itertools
import logging
import os
import shutil
import time
from typing import Dict

from django.conf import settings
from django.utils import timezone

import numpy as np
import pandas as pd
import psutil
from fancyimpute import IterativeSVD

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    CompendiumResult,
    CompendiumResultOrganismAssociation,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Organism,
    Pipeline,
    Sample,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import smashing_utils, utils

pd.set_option("mode.chained_assignment", None)


S3_COMPENDIA_BUCKET_NAME = get_env_variable("S3_COMPENDIA_BUCKET_NAME", "data-refinery")
BYTES_IN_GB = 1024 * 1024 * 1024
SMASHING_DIR = "/home/user/data_store/smashed/"
logger = get_and_configure_logger(__name__)
# DEBUG #
logger.setLevel(logging.getLevelName("DEBUG"))


def log_state(message, job_id, start_time=False):
    if logger.isEnabledFor(logging.DEBUG):
        process = psutil.Process(os.getpid())
        ram_in_GB = process.memory_info().rss / BYTES_IN_GB
        logger.debug(message, total_cpu=psutil.cpu_percent(), process_ram=ram_in_GB, job_id=job_id)

        if start_time:
            logger.debug("Duration: %s" % (time.time() - start_time), job_id=job_id)
        else:
            return time.time()


def _prepare_input(job_context: Dict) -> Dict:
    start_time = log_state("prepare input", job_context["job"].id)

    job_context["primary_organism"] = max(
        job_context["samples"], key=lambda organism: len(job_context["samples"][organism])
    )
    job_context["all_organisms"] = job_context["samples"].keys()
    all_samples = list(itertools.chain(*job_context["samples"].values()))
    job_context["samples"] = {job_context["primary_organism"]: all_samples}

    # We'll store here all sample accession codes that didn't make it into the compendia
    # with the reason why not.
    job_context["filtered_samples"] = {}

    # submitter processed accession codes that will be filtered out
    submitter_processed_samples = Sample.objects.filter(
        organism__name__in=job_context["all_organisms"],
        results__processor__name="Submitter-processed",
    ).values_list("accession_code", flat=True)

    job_context["submitter_processed_samples"] = set(submitter_processed_samples)

    job_context = smashing_utils.prepare_files(job_context)

    # Compendia jobs only run for one organism, so we know the only
    # key will be the organism name, unless of course we've already failed.
    if job_context["job"].success is not False:
        job_context["organism_name"] = job_context["group_by_keys"][0]

        # TEMPORARY for iterating on compendia more quickly. Rather
        # than downloading the data from S3 each run we're just gonna
        # use the same directory every job.
        job_context["old_work_dir"] = job_context["work_dir"]
        job_context["work_dir"] = SMASHING_DIR + job_context["organism_name"] + "/"
        if not os.path.exists(job_context["work_dir"]):
            os.makedirs(job_context["work_dir"])

    log_state("prepare input done", job_context["job"].id, start_time)
    return job_context


def _prepare_frames(job_context: Dict) -> Dict:
    start_prepare_frames = log_state("start _prepare_frames", job_context["job"].id)

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
            job_context = smashing_utils.process_frames_for_key(key, input_files, job_context)
            # if len(job_context['all_frames']) < 1:
            # TODO: Enable this check?
    except Exception:
        raise utils.ProcessorJobError(
            "Could not prepare frames for compendia.",
            success=False,
            dataset_id=job_context["dataset"].id,
            processor_job_id=job_context["job_id"],
            num_input_files=job_context["num_input_files"],
        )

    job_context["dataset"].success = True
    job_context["dataset"].save()

    log_state("end _prepare_frames", job_context["job"].id, start_prepare_frames)
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
    imputation_start = log_state("start perform imputation", job_context["job"].id)
    job_context["time_start"] = timezone.now()
    rnaseq_row_sums_start = log_state("start rnaseq row sums", job_context["job"].id)

    # We potentially can have a microarray-only compendia but not a RNASeq-only compendia
    log2_rnaseq_matrix = None
    if job_context["rnaseq_matrix"] is not None:
        # Drop any genes that are entirely NULL in the RNA-Seq matrix
        job_context["rnaseq_matrix"] = job_context["rnaseq_matrix"].dropna(
            axis="columns", how="all"
        )

        # Calculate the sum of the lengthScaledTPM values for each row
        # (gene) of the rnaseq_matrix (rnaseq_row_sums)
        rnaseq_row_sums = np.sum(job_context["rnaseq_matrix"], axis=1)

        log_state("end rnaseq row sums", job_context["job"].id, rnaseq_row_sums_start)
        rnaseq_decile_start = log_state("start rnaseq decile", job_context["job"].id)

        # Calculate the 10th percentile of rnaseq_row_sums
        rnaseq_tenth_percentile = np.percentile(rnaseq_row_sums, 10)

        log_state("end rnaseq decile", job_context["job"].id, rnaseq_decile_start)
        drop_start = log_state("drop all rows", job_context["job"].id)
        # Drop all rows in rnaseq_matrix with a row sum < 10th
        # percentile of rnaseq_row_sums; this is now
        # filtered_rnaseq_matrix
        # TODO: This is probably a better way to do this with `np.where`
        rows_to_filter = []
        for (x, sum_val) in rnaseq_row_sums.items():
            if sum_val < rnaseq_tenth_percentile:
                rows_to_filter.append(x)

        del rnaseq_row_sums

        log_state("actually calling drop()", job_context["job"].id)

        filtered_rnaseq_matrix = job_context.pop("rnaseq_matrix").drop(rows_to_filter)

        del rows_to_filter

        log_state("end drop all rows", job_context["job"].id, drop_start)
        log2_start = log_state("start log2", job_context["job"].id)

        # log2(x + 1) transform filtered_rnaseq_matrix; this is now log2_rnaseq_matrix
        filtered_rnaseq_matrix_plus_one = filtered_rnaseq_matrix + 1
        log2_rnaseq_matrix = np.log2(filtered_rnaseq_matrix_plus_one)
        del filtered_rnaseq_matrix_plus_one
        del filtered_rnaseq_matrix

        log_state("end log2", job_context["job"].id, log2_start)
        cache_start = log_state("start caching zeroes", job_context["job"].id)

        # Cache our RNA-Seq zero values
        cached_zeroes = {}
        for column in log2_rnaseq_matrix.columns:
            cached_zeroes[column] = log2_rnaseq_matrix.index[
                np.where(log2_rnaseq_matrix[column] == 0)
            ]

        # Set all zero values in log2_rnaseq_matrix to NA, but make sure
        # to keep track of where these zeroes are
        log2_rnaseq_matrix[log2_rnaseq_matrix == 0] = np.nan

        log_state("end caching zeroes", job_context["job"].id, cache_start)

    outer_merge_start = log_state("start outer merge", job_context["job"].id)

    # Perform a full outer join of microarray_matrix and
    # log2_rnaseq_matrix; combined_matrix
    if log2_rnaseq_matrix is not None:
        combined_matrix = job_context.pop("microarray_matrix").merge(
            log2_rnaseq_matrix, how="outer", left_index=True, right_index=True
        )
    else:
        logger.info("Building compendia with only microarray data.", job_id=job_context["job"].id)
        combined_matrix = job_context.pop("microarray_matrix")

    log_state("ran outer merge, now deleteing log2_rnaseq_matrix", job_context["job"].id)

    del log2_rnaseq_matrix

    log_state("end outer merge", job_context["job"].id, outer_merge_start)
    drop_na_genes_start = log_state("start drop NA genes", job_context["job"].id)

    # # Visualize Prefiltered
    # output_path = job_context['output_dir'] + "pre_filtered_" + str(time.time()) + ".png"
    # visualized_prefilter = visualize.visualize(combined_matrix.copy(), output_path)

    # Remove genes (rows) with <=70% present values in combined_matrix
    thresh = combined_matrix.shape[1] * 0.7  # (Rows, Columns)
    # Everything below `thresh` is dropped
    row_filtered_matrix = combined_matrix.dropna(axis="index", thresh=thresh)

    del combined_matrix
    del thresh

    log_state("end drop NA genes", job_context["job"].id, drop_na_genes_start)
    drop_na_samples_start = log_state("start drop NA samples", job_context["job"].id)

    # # Visualize Row Filtered
    # output_path = job_context['output_dir'] + "row_filtered_" + str(time.time()) + ".png"
    # visualized_rowfilter = visualize.visualize(row_filtered_matrix.copy(), output_path)

    # Remove samples (columns) with <50% present values in combined_matrix
    # XXX: Find better test data for this!
    col_thresh = row_filtered_matrix.shape[0] * 0.5
    row_col_filtered_matrix_samples = row_filtered_matrix.dropna(axis="columns", thresh=col_thresh)
    row_col_filtered_matrix_samples_index = row_col_filtered_matrix_samples.index
    row_col_filtered_matrix_samples_columns = row_col_filtered_matrix_samples.columns

    log_state("end drop NA genes", job_context["job"].id, drop_na_samples_start)
    replace_zeroes_start = log_state("start replace zeroes", job_context["job"].id)

    for sample_accession_code in row_filtered_matrix.columns:
        if sample_accession_code not in row_col_filtered_matrix_samples_columns:
            sample = Sample.objects.get(accession_code=sample_accession_code)
            sample_metadata = sample.to_metadata_dict()
            job_context["filtered_samples"][sample_accession_code] = {
                **sample_metadata,
                "reason": "Sample was dropped because it had less than 50% present values.",
                "experiment_accession_code": smashing_utils.get_experiment_accession(
                    sample.accession_code, job_context["dataset"].data
                ),
            }
        # submitter processed data has proven to be unreliable in the past
        # so we will filter it all out until we can check it one at a time
        # https://github.com/AlexsLemonade/refinebio/issues/2114
        if sample_accession_code in job_context["submitter_processed_samples"]:
            sample = Sample.objects.get(accession_code=sample_accession_code)
            sample_metadata = sample.to_metadata_dict()
            job_context["filtered_samples"][sample_accession_code] = {
                **sample_metadata,
                "reason": "Sample was dropped because it was submitter-processed.",
                "experiment_accession_code": smashing_utils.get_experiment_accession(
                    sample.accession_code, job_context["dataset"].data
                ),
            }
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
        except Exception:
            logger.warn("Error when replacing zero")
            continue

    log_state("end replace zeroes", job_context["job"].id, replace_zeroes_start)
    transposed_zeroes_start = log_state("start replacing transposed zeroes", job_context["job"].id)

    # Label our new replaced data
    combined_matrix_zero = row_col_filtered_matrix_samples
    del row_col_filtered_matrix_samples

    transposed_matrix_with_zeros = combined_matrix_zero.T
    del combined_matrix_zero

    # Remove -inf and inf
    # This should never happen, but make sure it doesn't!
    transposed_matrix = transposed_matrix_with_zeros.replace([np.inf, -np.inf], np.nan)
    del transposed_matrix_with_zeros

    log_state("end replacing transposed zeroes", job_context["job"].id, transposed_zeroes_start)

    # Store the absolute/percentages of imputed values
    matrix_sum = transposed_matrix.isnull().sum()
    percent = (matrix_sum / transposed_matrix.isnull().count()).sort_values(ascending=False)
    total_percent_imputed = sum(percent) / len(transposed_matrix.count())
    job_context["total_percent_imputed"] = total_percent_imputed
    logger.info("Total percentage of data to impute!", total_percent_imputed=total_percent_imputed)

    # Perform imputation of missing values with IterativeSVD (rank=10) on the
    # transposed_matrix; imputed_matrix
    svd_algorithm = job_context["dataset"].svd_algorithm
    if svd_algorithm != "NONE":
        svd_start = log_state("start SVD", job_context["job"].id)

        logger.info("IterativeSVD algorithm: %s" % svd_algorithm)
        svd_algorithm = str.lower(svd_algorithm)
        imputed_matrix = IterativeSVD(rank=10, svd_algorithm=svd_algorithm).fit_transform(
            transposed_matrix
        )

        svd_start = log_state("end SVD", job_context["job"].id, svd_start)
    else:
        imputed_matrix = transposed_matrix
        logger.info("Skipping IterativeSVD")
    del transposed_matrix

    untranspose_start = log_state("start untranspose", job_context["job"].id)

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
    job_context["organism"] = Organism.get_object_for_name(job_context["organism_name"])
    job_context["merged_no_qn"] = untransposed_imputed_matrix_df
    # output_path = job_context['output_dir'] + "compendia_no_qn_" + str(time.time()) + ".png"
    # visualized_merged_no_qn = visualize.visualize(untransposed_imputed_matrix_df.copy(),
    #                                               output_path)

    log_state("end untranspose", job_context["job"].id, untranspose_start)
    quantile_start = log_state("start quantile normalize", job_context["job"].id)

    # Perform the Quantile Normalization
    job_context = smashing_utils.quantile_normalize(job_context, ks_check=False)

    log_state("end quantile normalize", job_context["job"].id, quantile_start)

    # Visualize Final Compendia
    # output_path = job_context['output_dir'] + "compendia_with_qn_" + str(time.time()) + ".png"
    # visualized_merged_qn = visualize.visualize(job_context['merged_qn'].copy(), output_path)

    job_context["time_end"] = timezone.now()
    job_context["formatted_command"] = ["create_compendia.py"]
    log_state("end prepare imputation", job_context["job"].id, imputation_start)
    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    """
    Store and host the result as a ComputationalResult object.
    """
    result_start = log_state("start create result object", job_context["job"].id)
    result = ComputationalResult()
    result.commands.append(" ".join(job_context["formatted_command"]))
    result.is_ccdl = True
    # Temporary until we re-enable the QN test step.
    result.is_public = False
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    try:
        processor_key = "CREATE_COMPENDIA"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    # Write the compendia dataframe to a file
    job_context["csv_outfile"] = job_context["output_dir"] + job_context["organism_name"] + ".tsv"
    job_context["merged_qn"].to_csv(job_context["csv_outfile"], sep="\t", encoding="utf-8")

    organism_key = list(job_context["samples"].keys())[0]
    annotation = ComputationalResultAnnotation()
    annotation.result = result

    annotation.data = {
        "organism_id": job_context["samples"][organism_key][0].organism_id,
        "organism_name": job_context["organism_name"],
        "is_qn": False,
        "is_compendia": True,
        "samples": [sample.accession_code for sample in job_context["samples"][organism_key]],
        "num_samples": len(job_context["samples"][organism_key]),
        "experiment_accessions": [e.accession_code for e in job_context["experiments"]],
        "total_percent_imputed": job_context["total_percent_imputed"],
    }
    annotation.save()

    # Create the resulting archive
    final_zip_base = SMASHING_DIR + str(job_context["dataset"].pk) + "_compendia"
    # Copy LICENSE.txt and correct README.md files.
    if job_context["dataset"].quant_sf_only:
        readme_file = "/home/user/README_QUANT.md"
    else:
        readme_file = "/home/user/README_NORMALIZED.md"

    shutil.copy(readme_file, job_context["output_dir"] + "/README.md")
    shutil.copy("/home/user/LICENSE_DATASET.txt", job_context["output_dir"] + "/LICENSE.TXT")
    archive_path = shutil.make_archive(final_zip_base, "zip", job_context["output_dir"])

    archive_computed_file = ComputedFile()
    archive_computed_file.absolute_file_path = archive_path
    archive_computed_file.filename = archive_path.split("/")[-1]
    archive_computed_file.calculate_sha1()
    archive_computed_file.calculate_size()
    archive_computed_file.is_smashable = False
    archive_computed_file.is_qn_target = False
    archive_computed_file.result = result
    archive_computed_file.save()

    # Compendia Result Helpers
    primary_organism = Organism.get_object_for_name(job_context["primary_organism"])
    organisms = [
        Organism.get_object_for_name(organism) for organism in job_context["all_organisms"]
    ]
    compendium_version = (
        CompendiumResult.objects.filter(
            primary_organism=primary_organism, quant_sf_only=False
        ).count()
        + 1
    )
    # Save Compendia Result
    compendium_result = CompendiumResult()
    compendium_result.quant_sf_only = job_context["dataset"].quant_sf_only
    compendium_result.svd_algorithm = job_context["dataset"].svd_algorithm
    compendium_result.compendium_version = compendium_version
    compendium_result.result = result
    compendium_result.primary_organism = primary_organism
    compendium_result.save()

    # create relations to all organisms contained in the compendia

    compendium_result_organism_associations = []
    for compendium_organism in organisms:
        compendium_result_organism_association = CompendiumResultOrganismAssociation()
        compendium_result_organism_association.compendium_result = compendium_result
        compendium_result_organism_association.organism = compendium_organism
        compendium_result_organism_associations.append(compendium_result_organism_association)

    CompendiumResultOrganismAssociation.objects.bulk_create(compendium_result_organism_associations)

    job_context["compendium_result"] = compendium_result

    logger.info(
        "Compendium created!", archive_path=archive_path, organism_name=job_context["organism_name"]
    )

    # Upload the result to S3
    timestamp = str(int(time.time()))
    key = job_context["organism_name"] + "_" + str(compendium_version) + "_" + timestamp + ".zip"
    uploaded_to_s3 = archive_computed_file.sync_to_s3(S3_COMPENDIA_BUCKET_NAME, key)

    if not uploaded_to_s3:
        raise utils.ProcessorJobError(
            "Failed to upload compendia to S3",
            success=False,
            computed_file_id=archive_computed_file.id,
        )

    if settings.RUNNING_IN_CLOUD:
        archive_computed_file.delete_local_file()

    job_context["result"] = result
    job_context["success"] = True

    log_state("end create result object", job_context["job"].id, result_start)

    # TEMPORARY for iterating on compendia more quickly.
    # Reset this so the end_job does clean up the job's non-input-data stuff.
    job_context["work_dir"] = job_context["old_work_dir"]

    return job_context


def create_compendia(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.CREATE_COMPENDIA.value)
    job_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _prepare_input,
            _prepare_frames,
            _perform_imputation,
            smashing_utils.write_non_data_files,
            _create_result_objects,
            utils.end_job,
        ],
    )
    return job_context

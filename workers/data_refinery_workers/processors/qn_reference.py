import subprocess
import time
from typing import Dict

from django.utils import timezone

import pandas as pd

from data_refinery_common.enums import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Experiment,
    Pipeline,
)
from data_refinery_common.utils import queryset_iterator
from data_refinery_workers.processors import smashing_utils, utils

logger = get_and_configure_logger(__name__)


def _prepare_input(job_context: Dict) -> Dict:

    # We're going to use the smasher outside of the smasher.
    # I'm not crazy about this yet. Maybe refactor later,
    # but I need the data now.
    job_context = smashing_utils.prepare_files(job_context)

    # work_dir is already created by smasher._prepare_files
    outfile_base = job_context["work_dir"] + str(time.time()).split(".")[0]
    job_context["target_file"] = outfile_base + "_target.tsv"

    return job_context


def _build_qn_target(job_context: Dict) -> Dict:
    """Iteratively creates a QN target file, method described here:
    https://github.com/AlexsLemonade/refinebio/pull/1013
    """
    job_context["time_start"] = timezone.now()

    # Get the gene list from the first input
    (computed_file, _) = job_context["input_files"]["ALL"][0]
    computed_file_path = computed_file.get_synced_file_path()
    geneset_target_frame = smashing_utils._load_and_sanitize_file(computed_file_path)

    # Get the geneset
    geneset = set(geneset_target_frame.index.values)

    # Sequentially build the target
    sum_frame_input = {"index": list(geneset), "sum": [0 for value in geneset]}
    sum_frame = pd.DataFrame(data=sum_frame_input)
    sum_frame = sum_frame.set_index("index")

    # Read and sum all of the inputs
    num_valid_inputs = 0
    for file, sample in job_context["input_files"]["ALL"]:
        try:
            input_filepath = file.get_synced_file_path()
            input_frame = smashing_utils._load_and_sanitize_file(input_filepath)
        except Exception:
            logger.warn(
                "No file loaded for input file",
                exc_info=1,
                bad_file=file,
                num_valid_inputs_so_far=num_valid_inputs,
            )
            continue

        # If this input doesn't have the same geneset, we don't want it!
        if set(input_frame.index.values) != geneset:
            logger.warn(
                "Input frame doesn't match target geneset, skipping!",
                bad_file=file,
                target_geneset_len=len(geneset),
                bad_geneset_len=len(input_frame.index.values),
                geneset_difference=list(geneset ^ set(input_frame.index.values))[:3],
                num_valid_inputs_so_far=num_valid_inputs,
            )
            continue

        # We don't want to build QNs that have Nulls in them, so
        # filter out any samples that still have null genes at this
        # point.
        num_nulls = input_frame.isnull().sum(axis=0)[0]
        if num_nulls != 0:
            logger.warn(
                "Input frame contains NA values, skipping!",
                bad_file=file,
                number_of_NAs=num_nulls,
                num_valid_inputs_so_far=num_valid_inputs,
            )
            continue

        # Sort the input
        sample_name = list(input_frame.columns.values)[0]
        expressions = input_frame.sort_values(sample_name)

        # Add the sorted input to the sum frame
        sum_frame["sum"] = sum_frame["sum"] + expressions[sample_name].values

        # We'll divide by this later
        num_valid_inputs = num_valid_inputs + 1

    # Divide our summation by the number of inputs and save the resulting object and metadata
    sum_frame["sum"] = sum_frame["sum"] / num_valid_inputs
    job_context["time_end"] = timezone.now()
    job_context["sum_frame"] = sum_frame
    job_context["num_valid_inputs"] = num_valid_inputs
    job_context["geneset"] = list(geneset)

    # Write the file
    if job_context["create_results"]:
        sum_frame.to_csv(
            job_context["target_file"], index=False, header=False, sep="\t", encoding="utf-8"
        )
    job_context["formatted_command"] = "qn_reference.py"

    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    if not job_context["create_results"]:
        return job_context

    result = ComputationalResult()
    result.commands.append(" ".join(job_context["formatted_command"]))
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    try:
        processor_key = "QN_REFERENCE"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context["target_file"]
    computed_file.filename = job_context["target_file"].split("/")[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.is_smashable = False
    computed_file.is_qn_target = True
    computed_file.result = result
    computed_file.save()

    annotation = ComputationalResultAnnotation()
    annotation.result = result
    annotation.data = {
        "organism_id": job_context["samples"]["ALL"][0].organism_id,
        "is_qn": True,
        "platform_accession_code": job_context["samples"]["ALL"][0].platform_accession_code,
        "samples": [sample.accession_code for sample in job_context["samples"]["ALL"]],
        "geneset": str(job_context["geneset"]),
        "num_valid_inputs": job_context["num_valid_inputs"],
    }
    annotation.save()

    job_context["result"] = result
    job_context["computed_files"] = [computed_file]
    job_context["annotation"] = annotation
    job_context["success"] = True
    return job_context


def _update_caches(job_context: Dict) -> Dict:
    """Experiments have a cached value with the number of samples that have QN targets
    generated, this value should be updated after generating new QN targets."""
    if not job_context["create_results"]:
        # No new results were created, so there's no point in updating the cache
        return job_context

    organism = job_context["samples"]["ALL"][0].organism

    if job_context["result"]:
        # Associate the organism with the latest qn target that was generated for it
        organism.qn_target = job_context["result"]
        organism.save()

    unique_experiments = Experiment.objects.all().filter(organism=organism).distinct()

    for experiment in queryset_iterator(unique_experiments):
        experiment.update_num_samples()

    return job_context


def create_qn_reference(job_id: int, create_results=True) -> None:
    pipeline = Pipeline(name=PipelineEnum.QN_REFERENCE.value)
    job_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline, "create_results": create_results},
        [
            utils.start_job,
            _prepare_input,
            _build_qn_target,
            _create_result_objects,
            _update_caches,
            utils.end_job,
        ],
    )
    return job_context

import logging
import os
import shutil
import time
from typing import Dict, List

from django.conf import settings
from django.utils import timezone

import psutil

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    CompendiumResult,
    ComputationalResult,
    ComputedFile,
    Organism,
    Pipeline,
    Sample,
)
from data_refinery_common.utils import FileUtils, get_env_variable
from data_refinery_workers.processors import smashing_utils, utils

S3_COMPENDIA_BUCKET_NAME = get_env_variable("S3_COMPENDIA_BUCKET_NAME", "data-refinery")
SMASHING_DIR = "/home/user/data_store/smashed/"

logger = get_and_configure_logger(__name__)
logger.setLevel(logging.getLevelName("DEBUG"))


def create_quantpendia(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.CREATE_QUANTPENDIA.value)
    job_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _make_dirs,
            _download_files,
            _add_metadata,
            _make_archive,
            _create_result_objects,
            _remove_job_dir,
            utils.end_job,
        ],
    )
    return job_context


@utils.cache_keys(
    "time_start", "num_samples", "time_end", "formatted_command", work_dir_key="job_dir"
)
def _download_files(job_context: Dict) -> Dict:
    job_context["filtered_samples"] = {}
    job_context["time_start"] = timezone.now()

    num_samples = 0
    for key, samples in job_context["samples"].items():
        outfile_dir = job_context["output_dir"] + key + "/"
        os.makedirs(outfile_dir, exist_ok=True)

        logger.debug(
            "Downloading quant.sf files for quantpendia.",
            accession_code=key,
            job_id=job_context["job_id"],
            **get_process_stats()
        )

        # download quant.sf files directly into the dataset folder
        num_samples += smashing_utils.sync_quant_files(
            outfile_dir, samples, job_context["filtered_samples"]
        )

    job_context["num_samples"] = num_samples
    job_context["time_end"] = timezone.now()
    job_context["formatted_command"] = "create_quantpendia.py"

    logger.debug(
        "Finished downloading quant.sf files for quantpendia.",
        job_id=job_context["job_id"],
        total_downloaded_files=num_samples,
        **get_process_stats()
    )

    return job_context


@utils.cache_keys("metadata", work_dir_key="job_dir")
def _add_metadata(job_context: Dict) -> Dict:
    logger.debug(
        "Writing metadata for quantpendia.", job_id=job_context["job_id"], **get_process_stats()
    )
    smashing_utils.write_non_data_files(job_context)
    shutil.copy("/home/user/README_QUANT.md", job_context["output_dir"] + "/README.md")
    return job_context


@utils.cache_keys("archive_path", work_dir_key="job_dir")
def _make_archive(job_context: Dict):
    compendia_organism = _get_organisms(job_context["samples"]).first()
    final_zip_base = job_context["job_dir"] + compendia_organism.name + "_rnaseq_compendia"

    logger.debug(
        "Generating archive.",
        job_id=job_context["job_id"],
        organism_name=compendia_organism.name,
        **get_process_stats()
    )
    archive_path = shutil.make_archive(final_zip_base, "zip", job_context["output_dir"])
    logger.debug(
        "Quantpendia zip file generated.",
        job_id=job_context["job_id"],
        organism_name=compendia_organism.name,
        **get_process_stats()
    )

    return {**job_context, "archive_path": archive_path}


def _create_result_objects(job_context: Dict) -> Dict:
    """
    Store and host the result as a ComputationalResult object.
    """
    archive_path = job_context["archive_path"]
    compendia_organism = _get_organisms(job_context["samples"]).first()
    compendia_version = _get_next_compendia_version(compendia_organism)

    result = ComputationalResult()
    result.commands.append(" ".join(job_context["formatted_command"]))
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    try:
        processor_key = "CREATE_QUANTPENDIA"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    archive_computed_file = ComputedFile()
    archive_computed_file.absolute_file_path = archive_path
    archive_computed_file.filename = FileUtils.get_filename(archive_path)
    archive_computed_file.calculate_sha1()
    archive_computed_file.calculate_size()
    archive_computed_file.is_smashable = False
    archive_computed_file.is_qn_target = False
    archive_computed_file.result = result
    archive_computed_file.is_compendia = True
    archive_computed_file.quant_sf_only = True
    archive_computed_file.compendia_organism = compendia_organism
    archive_computed_file.compendia_version = compendia_version
    archive_computed_file.save()

    compendium_result = CompendiumResult()
    compendium_result.quant_sf_only = True
    compendium_result.result = result
    compendium_result.primary_organism = compendia_organism
    compendium_result.compendium_version = compendia_version
    compendium_result.save()

    logger.info(
        "Quantpendia created! Uploading to S3.",
        job_id=job_context["job_id"],
        archive_path=archive_path,
        organism_name=compendia_organism.name,
        **get_process_stats()
    )

    # Upload the result to S3
    timestamp = str(int(time.time()))
    s3_key = compendia_organism.name + "_" + str(compendia_version) + "_" + timestamp + ".zip"
    uploaded_to_s3 = archive_computed_file.sync_to_s3(S3_COMPENDIA_BUCKET_NAME, s3_key)

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

    return job_context


def _remove_job_dir(job_context: Dict):
    """ remove the directory when the job is successful. At this point
    the quantpendia was already zipped and uploaded. """
    # don't remove the files when running locally or for tests
    if settings.RUNNING_IN_CLOUD:
        shutil.rmtree(job_context["job_dir"], ignore_errors=True)
    return job_context


def _make_dirs(job_context: Dict):
    dataset_id = str(job_context["dataset"].pk)
    job_context["job_dir"] = "/home/user/data_store/smashed/" + dataset_id + "/"
    os.makedirs(job_context["job_dir"], exist_ok=True)
    job_context["output_dir"] = job_context["job_dir"] + "output/"
    os.makedirs(job_context["output_dir"], exist_ok=True)
    return job_context


def get_process_stats():
    BYTES_IN_GB = 1024 * 1024 * 1024
    process = psutil.Process(os.getpid())
    ram_in_GB = process.memory_info().rss / BYTES_IN_GB
    return {"total_cpu": psutil.cpu_percent(), "process_ram": ram_in_GB}


def _get_organisms(aggregated_samples: Dict[str, Sample]) -> List[Organism]:
    organisms = set()
    for key, samples in aggregated_samples.items():
        organism_ids = samples.values_list("organism__id", flat=True).distinct()
        organisms.update(organism_ids)

    return Organism.objects.filter(id__in=list(organisms))


def _get_next_compendia_version(organism: Organism) -> int:
    last_compendia = (
        ComputedFile.objects.filter(
            is_compendia=True, quant_sf_only=True, compendia_organism=organism
        )
        .order_by("-compendia_version")
        .first()
    )

    if last_compendia:
        return last_compendia.compendia_version + 1

    # otherwise this is the first compendia that we are generating
    return 1

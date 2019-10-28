import os
import logging
import shutil
import time
from django.utils import timezone
from django.conf import settings
from typing import Dict, List, Tuple
import psutil

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (ComputationalResult,
                                         ComputedFile,
                                         Organism,
                                         Pipeline,
                                         Sample)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import smashing_utils, utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
SMASHING_DIR = "/home/user/data_store/smashed/"

logger = get_and_configure_logger(__name__)
logger.setLevel(logging.getLevelName('DEBUG'))

def create_quantpendia(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.CREATE_QUANTPENDIA.value)
    job_context = utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                                     [utils.start_job,
                                      make_dirs,
                                      download_files,
                                      create_result_objects,
                                      remove_job_dir,
                                      utils.end_job])
    return job_context


def download_files(job_context: Dict) -> Dict:
    job_context['time_start'] = timezone.now()

    num_samples = 0
    for key, samples in job_context['samples'].items():
        outfile_dir = job_context['output_dir'] + key + '/'
        os.makedirs(outfile_dir, exist_ok=True)

        logger.debug("Downloading quant.sf files for quantpendia.",
                    accession_code=key,
                    job_id=job_context['job_id'],
                    **get_process_stats())

        # download quant.sf files directly into the dataset folder
        num_samples += smashing_utils.sync_quant_files(outfile_dir, samples)

    job_context['num_samples'] = num_samples
    job_context['time_end'] = timezone.now()
    job_context['formatted_command'] = "create_quantpendia.py"

    logger.debug("Finished downloading quant.sf files for quantpendia.",
                job_id=job_context['job_id'],
                total_downloaded_files=num_samples,
                **get_process_stats())

    return job_context


def create_result_objects(job_context: Dict) -> Dict:
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
        processor_key = "CREATE_QUANTPENDIA"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)
    result.save()

    compendia_organism = _get_organisms(job_context['samples']).first()

    # Create the resulting archive
    smashing_utils.write_non_data_files(job_context)
    final_zip_base = job_context['job_dir'] + compendia_organism.name + "_rnaseq_compendia"
    shutil.copy("/home/user/README_QUANT.md", job_context["output_dir"] + "/README.md")

    archive_path = shutil.make_archive(final_zip_base, 'zip', job_context["output_dir"])
    compendia_version = _get_next_compendia_version(compendia_organism)

    archive_computed_file = ComputedFile()

    archive_computed_file.absolute_file_path = archive_path
    archive_computed_file.filename = archive_path.split('/')[-1]
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

    logger.info("Quantpendia created!",
                archive_path=archive_path,
                organism_name=compendia_organism.name)

    # Upload the result to S3
    timestamp = str(int(time.time()))
    s3_key = compendia_organism.name + "_" + str(compendia_version) + "_" + timestamp + ".zip"
    archive_computed_file.sync_to_s3(S3_BUCKET_NAME, s3_key)

    job_context['result'] = result
    job_context['success'] = True

    return job_context


def remove_job_dir(job_context: Dict):
    """ remove the directory when the job is successful. At this point
    the quantpendia was already zipped and uploaded. """
    # don't remove the files when running locally or for tests
    if settings.RUNNING_IN_CLOUD:
        shutil.rmtree(job_context["job_dir"], ignore_errors=True)
    return job_context

def make_dirs(job_context: Dict):
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
    return { 'total_cpu': psutil.cpu_percent(), 'process_ram': ram_in_GB }


def _get_organisms(aggregated_samples: Dict[str, Sample]) -> List[Organism]:
    organisms = set()
    for key, samples in aggregated_samples.items():
        organism_ids = samples.values_list('organism__id', flat=True).distinct()
        organisms.update(organism_ids)

    return Organism.objects.filter(id__in=list(organisms))


def _get_next_compendia_version(organism: Organism) -> int:
    last_compendia = ComputedFile.objects\
        .filter(is_compendia=True, quant_sf_only=True, compendia_organism=organism)\
        .order_by('-compendia_version').first()

    if last_compendia:
        return last_compendia.compendia_version + 1

    # otherwise this is the first compendia that we are generating
    return 1

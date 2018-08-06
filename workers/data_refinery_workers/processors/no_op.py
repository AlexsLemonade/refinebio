import csv
import os
import shutil
import boto3

import subprocess
from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    SampleResultAssociation,
    SampleComputedFileAssociation,
    SampleAnnotation,
    Processor,
    Pipeline
)
from data_refinery_common.utils import get_env_variable, get_internal_microarray_accession
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """A processor which takes externally-processed sample data and makes it smashable.
    """
    try:
        original_file = job_context["original_files"][0]
        sample0 = job_context['samples'][0]
        if sample0.manufacturer == 'ILLUMINA':
            job_context["is_illumina"] = True
        else:
            job_context["is_illumina"] = False

        # Create the output directory and path
        job_context["input_file_path"] = original_file.absolute_file_path
        base_directory, file_name = original_file.absolute_file_path.rsplit('/', 1)
        os.makedirs(base_directory + '/processed/', exist_ok=True)
        job_context["output_file_path"] = base_directory + '/processed/' + file_name

        # Make sure header column is correct
        with open(job_context["input_file_path"], 'r', encoding='utf-8') as tsv_in:
            tsv_in = csv.reader(tsv_in, delimiter='\t')
            for line in tsv_in:
                row = line
                joined = ''.join(row)
                break

        # We want to make sure that all of our columns for conversion and smashing
        # use the same column name for gene identifiers for later lookup. ID_REF
        # is the most common, so we use that.
        if 'ID_REF' not in joined and not job_context["is_illumina"]:
            try:
                float(row[1])
                # Okay, there's no header so can just prepend to the file.
                with open(job_context["input_file_path"], 'r', encoding='utf-8') as f:
                    all_content = f.read()

                job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
                with open(job_context["input_file_path"], 'w+', encoding='utf-8') as f:
                    f.seek(0, 0)
                    f.write('ID_REF\tVALUE' + '\n' + all_content)

            except ValueError:
                # There is already a header row. Let's replace it with ID_Ref
                with open(job_context["input_file_path"], 'r', encoding='utf-8') as f:
                    all_content = f.read()

                job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
                with open(job_context["input_file_path"], 'w+', encoding='utf-8') as f:
                    f.seek(0, 0)
                    f.write('ID_REF\tVALUE' + '\n' + ("\n".join(all_content[1:])))
                job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
        else:
            if job_context["is_illumina"]:
                job_context['column_name'] = row[0]

        # Platform
        job_context["platform"] = job_context["samples"][0].platform_accession_code
        job_context["internal_accession"] = get_internal_microarray_accession(job_context["platform"])
        if not job_context["internal_accession"]:
            logger.error("Failed to find internal accession for code", code=job_context["platform"])
            job_context['success'] = False
    except Exception as e:
        logger.exception("Failed to prepare for NO_OP job.", job_id=job_context['job'].id)
        job_context["job"].failure_reason = str(e)
        job_context['success'] = False

    return job_context

def _convert_genes(job_context: Dict) -> Dict:
    """ Dispatches to the appropriate gene converter"""

    job_context["success"] = True

    sample0 = job_context['samples'][0]
    if sample0.manufacturer == 'ILLUMINA':
        return _convert_illumina_genes(job_context)
    else:
        return _convert_affy_genes(job_context)

def _convert_affy_genes(job_context: Dict) -> Dict:
    """ Convert to Ensembl genes if we can"""

    gene_index_path = "/home/user/gene_indexes/" + job_context["internal_accession"] + ".tsv.gz"
    if not os.path.exists(gene_index_path):
        logger.error("Missing gene index file for platform!",
            platform=job_context["internal_accession"],
            job_id=job_context["job_id"])
        job_context["failure_reason"] = "Missing gene index for " + job_context['internal_accession']
        job_context["success"] = False
        return job_context

    job_context['script_name'] = "gene_convert.R"
    try:
        result = subprocess.check_output([
                "/usr/bin/Rscript",
                "--vanilla",
                "/home/user/data_refinery_workers/processors/" + job_context['script_name'],
                "--platform", job_context["internal_accession"],
                "--inputFile", job_context['input_file_path'],
                "--outputFile", job_context['output_file_path']
            ])
    except Exception as e:
        error_template = ("Encountered error in R code while running gene_convert.R"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(job_context['input_file_path'], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        return job_context

    job_context["success"] = True
    return job_context


def _convert_illumina_genes(job_context: Dict) -> Dict:
    """ Convert to Ensembl genes if we can"""

    all_databases = {
        'HOMO_SAPIENS': [
            'illuminaHumanv1',
            'illuminaHumanv2',
            'illuminaHumanv3',
            'illuminaHumanv4',
        ],
        'MUS_MUSCULUS': [
            'illuminaMousev1',
            'illuminaMousev1p1',
            'illuminaMousev2',
        ],
        'RATTUS_NORVEGICUS': [
            'illuminaRatv1'
        ]
    }

    sample0 = job_context['samples'][0]
    databases = all_databases[sample0.organism.name]

    # Loop over all of the possible platforms and find the one with the best match.
    highest = 0.0
    high_mapped_percent = 0.0
    high_db = None
    for platform in databases:
        try:
            result = subprocess.check_output([
                    "/usr/bin/Rscript",
                    "--vanilla",
                    "/home/user/data_refinery_workers/processors/detect_database.R",
                    "--platform", platform,
                    "--inputFile", job_context['input_file_path'],
                    "--column", job_context.get('column_name', "Reporter Identifier")
                ])

            results = result.decode().split('\n')
            cleaned_result = float(results[0].strip())

            if cleaned_result > highest:
                highest = cleaned_result
                high_db = platform
                high_mapped_percent = float(results[1].strip())

        except Exception as e:
            logger.error("Could not detect database for file!",
                file=job_context['input_file_path'],
                platform=platform
            )
            logger.exception(e, context=job_context)
            continue

    # Record our sample detection outputs for every sample.
    for sample in job_context['samples']:
        sa = SampleAnnotation()
        sa.sample = sample
        sa.data = {
            "detected_platform": high_db,
            "detection_percentage": highest,
            "mapped_percentage": high_mapped_percent
        }
        sa.save()

    job_context['script_name'] = "gene_convert_illumina.R"
    try:
        result = subprocess.check_output([
                "/usr/bin/Rscript",
                "--vanilla",
                "/home/user/data_refinery_workers/processors/" + job_context['script_name'],
                "--platform", high_db,
                "--inputFile", job_context['input_file_path'],
                "--outputFile", job_context['output_file_path']
            ])
    except Exception as e:
        error_template = ("Encountered error in R code while running gene_convert_illumina.R"
                          " pipeline during processing of {0}: {1}")
        error_message = error_template.format(job_context['input_file_path'], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        return job_context

    job_context["success"] = True
    return job_context


def _create_result(job_context: Dict) -> Dict:
    """ Create the actual Result object"""

    # This is a NO-OP, but we make a ComputationalResult regardless.
    result = ComputationalResult()
    result.commands.append(job_context['script_name'])
    result.is_ccdl = True
    result.pipeline = "Submitter-processed"  # TODO: should be removed
    try:
        result.processor = utils.find_processor("SUBMITTER_PROCESSED")
    except Exception as e:
        err_str = "Failed to set processor: %s" % e
        logger.error(err_str)
        job_context["job"].failure_reason = err_str
        job_context["success"] = False
        return job_context

    result.save()
    job_context['pipeline'].steps.append(result.id)

    # Create a ComputedFile for the original file,
    # sync it S3 and save it.
    try:
        computed_file = ComputedFile()
        computed_file.absolute_file_path = job_context["output_file_path"]
        computed_file.filename = job_context['output_file_path'].split('/')[-1]
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        # computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
        # TODO here: delete local file after S3 sync
        computed_file.save()
    except Exception:
        logger.error("Exception caught while moving file %s",
                         raw_path,
                         processor_job=job_context["job_id"])

        failure_reason = "Exception caught while moving file {}".format(file.name)
        job_context["job"].failure_reason = failure_reason
        job_context["success"] = False
        return job_context

    for sample in job_context['samples']:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample,
            computed_file=computed_file)

    logger.info("Created %s", result)
    job_context["success"] = True
    return job_context

def no_op_processor(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.NO_OP.value)
    return utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                              [utils.start_job,
                               _prepare_files,
                               _convert_genes,
                               _create_result,
                               utils.end_job])

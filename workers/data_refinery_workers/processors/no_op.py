import boto3
import csv
import os
import shutil
import subprocess

from typing import Dict

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Pipeline,
    Processor,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable, get_internal_microarray_accession
from data_refinery_workers.processors import utils


logger = get_and_configure_logger(__name__)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")

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

        # All files for the job are in the same directory.
        job_context["work_dir"] = LOCAL_ROOT_DIR + "/" + "processor_job_" + str(job_context["job_id"]) + "/"
        try:
            os.makedirs(job_context["work_dir"])
        except Exception as e:
            logger.exception("Could not create work directory for processor job.",
                             job_context=job_context)
            job_context["job"].failure_reason = str(e)
            job_context["success"] = False
            return job_context

        # Create the output directory and path
        job_context["input_file_path"] = original_file.absolute_file_path
        new_filename = "gene_converted_" + original_file.filename
        job_context["output_file_path"] = job_context["work_dir"] + new_filename

        try:
            # Make sure header row is correct
            with open(job_context["input_file_path"], 'r', encoding='utf-8') as tsv_in:
                tsv_in = csv.reader(tsv_in, delimiter='\t')
                for line in tsv_in:
                    row = line
                    joined = ''.join(row)
                    break
        except Exception as e:
            logger.exception("Expected text file, was this in fact a binary file?",
                             job_context=job_context,
                             file=job_context["input_file_path"])

            job_context["job"].failure_reason = str(e)
            job_context["job"].no_retry = True
            job_context["success"] = False
            return job_context

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
            except Exception as e:
                logger.exception("Unable to read input file or header row.",
                    input_file_path=job_context["input_file_path"])
                job_context['job'].failure_reason = str(e)
                job_context['success'] = False
                job_context["job"].no_retry = True
                return job_context
        else:
            if job_context["is_illumina"]:
                try:
                    float(row[1])
                    # If we're here, there's no header. We're gonna have to fudge it a little bit.

                    # These exist, but they're too hard for us to handle right now.
                    if len(row) > 3:
                        e_msg = "We found an Illumina file to NO_OP that we're not set up to process yet. Tell Rich!"
                        logger.error(e_msg,
                            job_id = job_context["job"].pk,
                            file=job_context["input_file_path"]
                        )
                        job_context['success'] = False
                        job_context["job"].failure_reason = str(e_msg)
                        job_context["job"].no_retry = True
                        return job_context

                    # Okay, there's no header so can just prepend to the file.
                    with open(job_context["input_file_path"], 'r', encoding='utf-8') as f:
                        all_content = f.read()
                    job_context["input_file_path"] = job_context["input_file_path"] + ".fixed"
                    with open(job_context["input_file_path"], 'w+', encoding='utf-8') as f:
                        f.seek(0, 0)
                        if len(row) == 2:
                            f.write('Reporter Identifier\tVALUE' + '\n' + all_content)
                        elif len(row) == 3:
                            f.write('Reporter Identifier\tVALUE\tDetection Pval' + '\n' + all_content)

                    job_context['column_name'] = 'Reporter Identifier'
                except ValueError:
                    # Okay, there's a header column, we're good.
                    job_context['column_name'] = row[0]

        # Platform
        job_context["platform"] = job_context["samples"][0].platform_accession_code
        job_context["internal_accession"] = get_internal_microarray_accession(job_context["platform"])
        if not job_context["internal_accession"]:
            logger.error("Failed to find internal accession for code", code=job_context["platform"])
            job_context["job"].failure_reason = "Failed to find internal accession for code" + job_context["platform"]
            job_context['success'] = False
            job_context["job"].no_retry = True
    except Exception as e:
        logger.exception("Failed to prepare for NO_OP job.", job_id=job_context['job'].id)
        job_context["job"].failure_reason = str(e)
        job_context['success'] = False
        job_context["job"].no_retry = True

    return job_context

def _convert_genes(job_context: Dict) -> Dict:
    """ Dispatches to the appropriate gene converter"""

    sample0 = job_context['samples'][0]
    if sample0.manufacturer == 'ILLUMINA':
        return _convert_illumina_genes(job_context)
    else:
        return _convert_affy_genes(job_context)

def _convert_affy_genes(job_context: Dict) -> Dict:
    """ Convert to Ensembl genes if we can"""

    if 'internal_accession' not in job_context or not job_context['internal_accession']:
        error_msg = "Told to convert AFFY genes without an internal_accession - how did this happen?"
        logger.error(error_msg, job_context=job_context)
        job_context["job"].failure_reason = str(error_msg)
        job_context['success'] = False
        job_context["job"].no_retry = True
        return job_context

    gene_index_path = "/home/user/gene_indexes/" + job_context["internal_accession"] + ".tsv.gz"
    if not os.path.exists(gene_index_path):
        logger.error("Missing gene index file for platform!",
            platform=job_context["internal_accession"],
            job_id=job_context["job_id"])
        job_context["job"].failure_reason = "Missing gene index for " + job_context['internal_accession']
        job_context["job"].no_retry = True
        job_context["success"] = False
        return job_context

    job_context['script_name'] = "gene_convert.R"
    try:
        result = subprocess.check_output([
            "/usr/bin/Rscript",
            "--vanilla", # Shut up Rscript
            "/home/user/data_refinery_workers/processors/" + job_context['script_name'],
            "--geneIndexPath", gene_index_path,
            "--inputFile", job_context['input_file_path'],
            "--outputFile", job_context['output_file_path']
        ], stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        error_template = "Status code {0} from {1}: {2}"
        error_message = error_template.format(e.returncode, job_context['script_name'], e.stderr)
        logger.error(error_message, job_context=job_context)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        job_context["job"].no_retry = True
        return job_context
    except Exception as e:
        error_template = ("Encountered error in R code while running {0}"
                          " pipeline during processing of {1}: {2}")
        error_message = error_template.format(job_context['script_name'],
            job_context['input_file_path'], str(e))
        logger.error(error_message, job_context=job_context)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        job_context["job"].no_retry = True
        return job_context

    # Quality control!
    # Related: https://github.com/AlexsLemonade/refinebio/issues/614
    # Related: GSM102671
    is_good_quality = check_output_quality(job_context['output_file_path'])
    if not is_good_quality:
        job_context["success"] = False
        job_context["job"].failure_reason = "NO_OP output failed quality control check. (" + job_context['output_file_path'] + ")"
        job_context["job"].no_retry = True
        return job_context
    else:
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
            logger.exception("Could not detect database for file!",
                             platform=platform,
                             job_context=job_context)
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
            ], stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        error_template = "Status code {0} from {1}: {2}"
        error_message = error_template.format(e.returncode, job_context['script_name'], e.stderr)
        logger.error(error_message, job_context=job_context)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        job_context["job"].no_retry = True
        return job_context
    except Exception as e:
        error_template = ("Encountered error in R code while running {0}"
                          " pipeline during processing of {1}: {2}")
        error_message = error_template.format(job_context['script_name'],
            job_context['input_file_path'], str(e))
        logger.error(error_message, job_context=job_context)
        job_context["job"].failure_reason = error_message
        job_context["success"] = False
        job_context["job"].no_retry = True
        return job_context

    job_context["success"] = True
    return job_context


def _create_result(job_context: Dict) -> Dict:
    """ Create the actual Result object"""

    # This is a NO-OP, but we make a ComputationalResult regardless.
    result = ComputationalResult()
    result.commands.append(job_context['script_name'])
    result.is_ccdl = True
    try:
        processor_key = "SUBMITTER_PROCESSED"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)

    result.save()
    job_context['pipeline'].steps.append(result.id)

    # Create a ComputedFile for the computed file,
    # sync it S3 and save it.
    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context["output_file_path"]
    computed_file.filename = job_context['output_file_path'].split('/')[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.result = result
    computed_file.is_smashable = True
    computed_file.is_qc = False
    computed_file.save()

    # utils.end_job will sync this to S3 for us.
    job_context["computed_files"] = [computed_file]

    for sample in job_context['samples']:
        assoc = SampleResultAssociation()
        assoc.sample = sample
        assoc.result = result
        assoc.save()

        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample,
            computed_file=computed_file)

    logger.debug("Created %s", result)
    job_context["success"] = True
    return job_context


def check_output_quality(output_file_path: str):
    """ Verify our output file meets spec """
    try:
        with open(output_file_path, 'r') as tsv_in:
            tsv_in = csv.reader(tsv_in, delimiter='\t')
            for i, row in enumerate(tsv_in):
                # We could change this to a header column checker
                if i == 0:
                    continue

                # If there are more than 2 columns,
                # we consider this bad data. (Likely old machine/processing!)
                if len(row) > 2:
                    return False
    except Exception as e:
        return False

    return True


def no_op_processor(job_id: int) -> None:
    pipeline = Pipeline(name=utils.PipelineEnum.NO_OP.value)
    return utils.run_pipeline({"job_id": job_id, "pipeline": pipeline},
                              [utils.start_job,
                               _prepare_files,
                               _convert_genes,
                               _create_result,
                               utils.end_job])

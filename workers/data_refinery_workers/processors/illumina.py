import csv
import multiprocessing
import os
import string
import subprocess
import warnings
from typing import Dict

import numpy as np
import pandas as pd
from django.utils import timezone

from data_refinery_common.job_lookup import PipelineEnum, ProcessorPipeline
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.message_queue import send_job
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    OriginalFile,
    OriginalFileSampleAssociation,
    Pipeline,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
logger = get_and_configure_logger(__name__)


def _prepare_files(job_context: Dict) -> Dict:
    """Adds the keys "input_file_path" and "output_file_path" to
    job_context so everything is prepared for processing.
    """
    # All files for the job are in the same directory.
    job_context["work_dir"] = (
        LOCAL_ROOT_DIR + "/" + "processor_job_" + str(job_context["job_id"]) + "/"
    )
    os.makedirs(job_context["work_dir"], exist_ok=True)

    original_file = job_context["original_files"][0]
    sanitized_filename = original_file.absolute_file_path.split("/")[-1] + ".sanitized"
    job_context["input_file_path"] = job_context["work_dir"] + sanitized_filename

    new_filename = original_file.absolute_file_path.split("/")[-1].replace(".txt", ".PCL")
    job_context["output_file_path"] = job_context["work_dir"] + new_filename

    # Sanitize this file so R doesn't choke.
    # Some have comments, some have non-comment-comments.
    with open(original_file.absolute_file_path, "r") as file_input:
        with open(job_context["input_file_path"], "w") as file_output:
            for line in file_input:
                if (
                    "#" not in line
                    and line.strip() != ""
                    and line != "\n"
                    and "\t" in line
                    and line[0:3].upper() != "GSM"
                    and line[0] != "'"
                    and line[0] != '"'
                    and line[0] != "!"
                    and line[0] != "/"
                    and line[0] != "<"
                    and line[0] != "\t"
                ):
                    file_output.write(line)

    return job_context


def _detect_columns(job_context: Dict) -> Dict:
    """ Detect which columns match to which inputs.

    Related: https://github.com/AlexsLemonade/refinebio/issues/86#issuecomment-379308817

    We need to find:

        First column should be ID_REF or PROBE_ID and the type should be string.
        Detection Pval column
        Expression column (contains sample title and NOT 'BEAD')

        Header examples:
            ['ID_REF', 'LV-C&si-Control-1', 'Detection Pval',
            'LV-C&si-Control-2', 'Detection Pval', 'LV-C&si-Control-3', 'Detection
            Pval', 'LV-C&si-EZH2-1', 'Detection Pval', 'LV-C&si-EZH2-2', 'Detection
            Pval', 'LV-C&si-EZH2-3', 'Detection Pval', 'LV-EZH2&si-EZH2-1',
            'Detection Pval', 'LV-EZH2&si-EZH2-2', 'Detection Pval', 'LV-EZH2&si-
            EZH2-3', 'Detection Pval', 'LV-T350A&si-EZH2-1', 'Detection Pval', 'LV-
            T350A&si-EZH2-2', 'Detection Pval', 'LV-T350A&si-EZH2-3', 'Detection
            Pval']

    Adds the following keys to job_context:
        columnIds: the identifiers of columns which contain expression data
        probeId: which is the value of the column containing the probe identifiers.
        detectionPval: a string which identifies Pvalue columns
    """
    try:
        input_file = job_context["input_file_path"]
        headers = None
        with open(input_file, "r") as tsv_in:
            tsv_in = csv.reader(tsv_in, delimiter="\t")
            for row in tsv_in:
                headers = row
                break

        # Ex GSE45331_non-normalized.txt
        predicted_header = 0
        if headers[0].upper() in ["TARGETID", "TARGET_ID"]:
            predicted_header = 1

        # First the probe ID column
        if headers[predicted_header].upper() not in [
            "ID_REF",
            "PROBE_ID",
            "IDREF",
            "PROBEID",
            "REF_ID",
            "REFID",
            "IDPROBE",
            "ID_PROBE",
        ]:
            job_context["job"].failure_reason = (
                "Could not find any ID column in headers "
                + str(headers)
                + " for file "
                + job_context["input_file_path"]
            )
            job_context["success"] = False
            return job_context
        else:
            job_context["probeId"] = headers[predicted_header]

        # Then the detection Pvalue string, which is always(?) some form of 'Detection Pval'
        for header in headers:
            if "DETECTION_PVAL" in header.upper().replace(" ", "_"):
                # Could be <SAMPLE_ID>.Detection Pval, etc
                if "." in header:
                    job_context["detectionPval"] = header.split(".")[-1]
                else:
                    job_context["detectionPval"] = header
                break
        else:
            job_context["job"].failure_reason = "Could not detect PValue column!"
            job_context["success"] = False
            job_context["job"].no_retry = True
            return job_context

        # Then, finally, create an absolutely bonkers regular expression
        # which will explicitly hit on any sample which contains a sample
        # ID _and_ ignores the magical word 'BEAD', etc. Great!
        column_ids = ""
        for sample in job_context["samples"]:
            for offset, header in enumerate(headers, start=1):

                if sample.title == header:
                    column_ids = column_ids + str(offset) + ","
                    continue

                # Sometimes the title might actually be in the description field.
                # To find this, look in all the related SampleAnnotations.
                # Since there are multiple annotations, we need to break early before continuing.
                # Related: https://github.com/AlexsLemonade/refinebio/issues/499
                continue_me = False
                for annotation in sample.sampleannotation_set.filter(is_ccdl=False):
                    try:
                        if annotation.data.get("description", "")[0] == header:
                            column_ids = column_ids + str(offset) + ","
                            continue_me = True
                            break
                    except Exception:
                        pass
                if continue_me:
                    # Treat the header as the real title, as we will need it later.
                    sample.title = header
                    sample.save()
                    continue

                if header.upper().replace(" ", "_") == "RAW_VALUE":
                    column_ids = column_ids + str(offset) + ","
                    continue

                if (
                    sample.title.upper() in header.upper()
                    and "BEAD" not in header.upper()
                    and "NARRAYS" not in header.upper()
                    and "ARRAY_STDEV" not in header.upper()
                    and "PVAL" not in header.upper().replace(" ", "").replace("_", "")
                ):
                    column_ids = column_ids + str(offset) + ","
                    continue

        for offset, header in enumerate(headers, start=1):
            if "AVG_Signal" in header:
                column_ids = column_ids + str(offset) + ","
                continue

        # Remove the trailing comma
        column_ids = column_ids[:-1]
        job_context["columnIds"] = column_ids
    except Exception as e:
        job_context["job"].failure_reason = str(e)
        job_context["success"] = False
        logger.exception(
            "Failed to extract columns in " + job_context["input_file_path"], exception=str(e)
        )
        job_context["job"].no_retry = True
        return job_context

    return job_context


def _detect_platform(job_context: Dict) -> Dict:
    """
    Determine the platform/database to process this sample with.
    They often provide something like "V2" or "V 2", but we don't trust them so we detect it ourselves.

    Related: https://github.com/AlexsLemonade/refinebio/issues/232
    """

    all_databases = {
        "HOMO_SAPIENS": [
            "illuminaHumanv1",
            "illuminaHumanv2",
            "illuminaHumanv3",
            "illuminaHumanv4",
        ],
        "MUS_MUSCULUS": ["illuminaMousev1", "illuminaMousev1p1", "illuminaMousev2",],
        "RATTUS_NORVEGICUS": ["illuminaRatv1"],
    }

    sample0 = job_context["samples"][0]
    databases = all_databases[sample0.organism.name]

    # Loop over all of the possible platforms and find the one with the best match.
    highest = 0.0
    high_mapped_percent = 0.0
    high_db = None
    for platform in databases:
        try:
            result = subprocess.check_output(
                [
                    "/usr/bin/Rscript",
                    "--vanilla",
                    "/home/user/data_refinery_workers/processors/detect_database.R",
                    "--platform",
                    platform,
                    "--inputFile",
                    job_context["input_file_path"],
                    "--column",
                    job_context["probeId"],
                ]
            )

            results = result.decode().split("\n")
            cleaned_result = float(results[0].strip())

            if cleaned_result > highest:
                highest = cleaned_result
                high_db = platform
                high_mapped_percent = float(results[1].strip())

        except Exception as e:
            logger.exception(e, processor_job_id=job_context["job"].id)
            continue

    # Record our sample detection outputs for every sample.
    for sample in job_context["samples"]:
        sa = SampleAnnotation()
        sa.sample = sample
        sa.is_ccdl = True
        sa.data = {
            "detected_platform": high_db,
            "detection_percentage": highest,
            "mapped_percentage": high_mapped_percent,
        }
        sa.save()

    # If the match is over 75%, record this and process it on that platform.
    if high_mapped_percent > 75.0:
        job_context["platform"] = high_db
    # The match percentage is too low - send this to the no-opper instead.
    else:

        logger.info("Match percentage too low, NO_OP'ing and aborting.", job=job_context["job_id"])

        processor_job = ProcessorJob()
        processor_job.pipeline_applied = "NO_OP"
        processor_job.volume_index = job_context["job"].volume_index
        processor_job.ram_amount = job_context["job"].ram_amount
        processor_job.save()

        assoc = ProcessorJobOriginalFileAssociation()
        assoc.original_file = job_context["original_files"][0]
        assoc.processor_job = processor_job
        assoc.save()

        try:
            send_job(ProcessorPipeline.NO_OP, processor_job)
        except Exception as e:
            # Nomad dispatch error, likely during local test.
            logger.error(e, job=processor_job)

        job_context["abort"] = True

    return job_context


def _run_illumina(job_context: Dict) -> Dict:
    """Processes an input TXT file to an output PCL file using a custom R script.

    Expects a job_context which has been pre-populated with inputs, outputs
    and the column identifiers which the R script needs for processing.
    """
    input_file_path = job_context["input_file_path"]

    try:
        job_context["time_start"] = timezone.now()

        formatted_command = [
            "/usr/bin/Rscript",
            "--vanilla",
            "/home/user/data_refinery_workers/processors/illumina.R",
            "--probeId",
            job_context["probeId"],
            "--expression",
            job_context["columnIds"],
            "--detection",
            job_context["detectionPval"],
            "--platform",
            job_context["platform"],
            "--inputFile",
            job_context["input_file_path"],
            "--outputFile",
            job_context["output_file_path"],
            "--cores",
            str(multiprocessing.cpu_count()),
        ]

        subprocess.check_output(formatted_command)

        job_context["formatted_command"] = " ".join(formatted_command)

        job_context["time_end"] = timezone.now()

    except Exception as e:
        error_template = (
            "Encountered error in R code while running illumina.R"
            " pipeline during processing of {0}: {1}"
        )
        error_message = error_template.format(job_context["input_file_path"], str(e))
        logger.error(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

    return job_context


def _create_result_objects(job_context: Dict) -> Dict:
    """ Create the ComputationalResult objects after a Scan run is complete """

    result = ComputationalResult()
    result.commands.append(job_context["formatted_command"])
    result.is_ccdl = True
    result.is_public = True
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    try:
        processor_key = "ILLUMINA_SCAN"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)

    result.save()
    job_context["pipeline"].steps.append(result.id)

    # Split the result into smashable subfiles
    big_tsv = job_context["output_file_path"]
    data = pd.read_csv(big_tsv, sep="\t", header=0, index_col=0)
    individual_files = []
    frames = np.split(data, len(data.columns), axis=1)
    for frame in frames:
        filename = (
            frame.columns.values[0].replace("&", "").replace("*", "").replace(";", "") + ".tsv"
        )
        frame_path = job_context["work_dir"] + filename
        frame.to_csv(frame_path, sep="\t", encoding="utf-8")

        # This needs to be the same as the ones in the job context!
        try:
            sample = job_context["samples"].get(title=frame.columns.values[0])
        except Sample.DoesNotExist:
            logger.error(
                "Could not find sample for column while splitting Illumina file.",
                title=frame.columns.values[0],
                processor_job=job_context["job_id"],
                file_path=big_tsv,
            )
            continue

        computed_file = ComputedFile()
        computed_file.absolute_file_path = frame_path
        computed_file.filename = frame_path.split("/")[-1]
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        computed_file.is_public = True
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.save()
        job_context["computed_files"].append(computed_file)

        SampleResultAssociation.objects.get_or_create(sample=sample, result=result)

        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample, computed_file=computed_file
        )

        individual_files.append(computed_file)

    logger.debug("Created %s", result)
    job_context["success"] = True
    job_context["individual_files"] = individual_files
    job_context["result"] = result

    return job_context


def illumina_to_pcl(job_id: int) -> None:
    pipeline = Pipeline(name=PipelineEnum.ILLUMINA.value)
    return utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _prepare_files,
            _detect_columns,
            _detect_platform,
            _run_illumina,
            _create_result_objects,
            utils.end_job,
        ],
    )

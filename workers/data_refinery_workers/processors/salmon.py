import json
import os
import re
import shutil
import subprocess
import tarfile
from typing import Dict, List

from django.conf import settings
from django.db import transaction
from django.utils import timezone

import boto3
import numpy as np
import pandas as pd
import untangle
from botocore.client import Config

from data_refinery_common.enums import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Experiment,
    OrganismIndex,
    Pipeline,
    Sample,
    SampleComputedFileAssociation,
    SampleResultAssociation,
)
from data_refinery_common.rna_seq import get_tximport_inputs_if_eligible
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import utils

# We have to set the signature_version to v4 since us-east-1 buckets require
# v4 authentication.
S3 = boto3.client("s3", config=Config(signature_version="s3v4"))
logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "processor_job_"
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")


def _set_job_prefix(job_context: Dict) -> Dict:
    """Sets the `job_dir_prefix` value in the job context object."""
    job_context["job_dir_prefix"] = JOB_DIR_PREFIX + str(job_context["job_id"])
    return job_context


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the file(s) from the raw directory to the temp directory.

    Also adds the keys "input_file_path" and "output_directory" to
    job_context so everything is prepared for processing. If the reads
    are paired then there will also be an "input_file_path_2" key
    added to job_context for the second read.
    """
    logger.debug("Preparing files..")

    # Create a directory specific to this processor job combo.
    # (A single sample could belong to multiple experiments, meaning
    # that it could be run more than once, potentially even at the
    # same time.)
    job_context["work_dir"] = os.path.join(LOCAL_ROOT_DIR, job_context["job_dir_prefix"]) + "/"
    os.makedirs(job_context["work_dir"], exist_ok=True)

    original_files = job_context["original_files"]
    job_context["input_file_path"] = original_files[0].absolute_file_path

    if not os.path.exists(job_context["input_file_path"]):
        logger.error(
            "Was told to process a non-existent file - why did this happen?",
            input_file_path=job_context["input_file_path"],
            processor_job=job_context["job_id"],
        )
        job_context["job"].failure_reason = "Missing input file: " + str(
            job_context["input_file_path"]
        )
        job_context["success"] = False
        return job_context

    if len(original_files) == 2:
        job_context["input_file_path_2"] = original_files[1].absolute_file_path
        if not os.path.exists(job_context["input_file_path_2"]):
            logger.error(
                "Was told to process a non-existent file2 - why did this happen?",
                input_file_path=job_context["input_file_path_2"],
                processor_job=job_context["job_id"],
            )
            job_context["job"].failure_reason = "Missing input file2: " + str(
                job_context["input_file_path_2"]
            )
            job_context["success"] = False
            return job_context

    # There should only ever be one per Salmon run
    sample = job_context["original_files"][0].samples.first()

    # This check was added to ensure that we don't process any RNA-Seq
    # samples from GEO, but for the time being we really don't want to
    # run salmon on anything that's not from SRA. See
    # https://github.com/AlexsLemonade/refinebio/issues/966 for more
    # information.
    if sample.technology != "RNA-SEQ" or sample.source_database != "SRA":
        failure_reason = (
            "The sample for this job either was not RNA-Seq or was not from the " "SRA database."
        )
        job_context["failure_reason"] = failure_reason
        logger.error(failure_reason, sample=sample, processor_job=job_context["job_id"])

        # We can't process it, we should mark it as such.
        sample.is_unable_to_be_processed = True
        sample.save()

        # No need to retry and fail more than once for this reason.
        job_context["success"] = False
        job_context["job"].failure_reason = failure_reason
        job_context["job"].no_retry = True
        return job_context

    # Detect that this is an SRA file from the source URL
    if ("ncbi.nlm.nih.gov" in job_context["original_files"][0].source_url) or (
        job_context["input_file_path"][-4:].upper() == ".SRA"
    ):
        new_input_file_path = os.path.join(job_context["work_dir"], original_files[0].filename)
        shutil.copyfile(job_context["input_file_path"], new_input_file_path)
        job_context["input_file_path"] = new_input_file_path
        job_context["sra_input_file_path"] = new_input_file_path

    if job_context.get("input_file_path_2", False):
        new_input_file_path = os.path.join(job_context["work_dir"], original_files[1].filename)
        shutil.copyfile(job_context["input_file_path_2"], new_input_file_path)
        job_context["input_file_path_2"] = new_input_file_path

    job_context["sample_accession_code"] = sample.accession_code
    job_context["sample"] = sample
    job_context["samples"] = []  # This will only be populated in the `tximport` job
    job_context["organism"] = job_context["sample"].organism
    job_context["success"] = True

    job_context["output_directory"] = job_context["work_dir"] + sample.accession_code + "_output/"
    os.makedirs(job_context["output_directory"], exist_ok=True)

    job_context["salmontools_directory"] = job_context["work_dir"] + "salmontools/"
    os.makedirs(job_context["salmontools_directory"], exist_ok=True)
    job_context["salmontools_archive"] = job_context["work_dir"] + "salmontools-result.tar.gz"

    timestamp = str(timezone.now().timestamp()).split(".")[0]
    job_context["output_archive"] = job_context["work_dir"] + "result-" + timestamp + ".tar.gz"

    job_context["computed_files"] = []
    job_context["smashable_files"] = []

    return job_context


def _determine_index_length_sra(job_context: Dict) -> Dict:
    """
    Use the sra-stat tool to determine length
    ex:
        sra-stat -x --statistics ERR1562482.sra
    """

    command_str = "sra-stat -x --statistics {sra_file}"
    formatted_command = command_str.format(sra_file=job_context["sra_input_file_path"])

    logger.debug(
        "Running sra-stat using the following shell command: %s",
        formatted_command,
        processor_job=job_context["job_id"],
    )

    completed_command = subprocess.run(
        formatted_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    respo = completed_command.stdout.decode().strip()

    try:
        stats = untangle.parse(respo)
    except ValueError:
        logger.error(
            "Unable to parse sra-stat output!", respo=str(respo), command=formatted_command
        )

    # Different SRA files can create different output formats, somehow.
    # This mess tries every output method we can to parse these stats.
    # If it's so messed up we don't know, default to short.
    try:
        bases_count = int(stats.Run.Bases["count"])
        reads_count = int(stats.Run.Statistics["nspots"]) * int(stats.Run.Statistics["nreads"])
        num_reads = int(stats.Run.Statistics["nreads"])
        job_context["index_length_raw"] = int(bases_count / reads_count)
    except Exception:
        try:
            num_reads = int(stats.Run.Statistics["nreads"])
            spot_count_mates = int(stats.Run["spot_count_mates"])
            base_count_bio_mates = int(stats.Run["base_count_bio_mates"])
            reads_count = spot_count_mates * int(stats.Run.Statistics["nreads"])
            job_context["index_length_raw"] = int(base_count_bio_mates / reads_count)
        except Exception:
            try:
                job_context["index_length_raw"] = int(stats.Run.Statistics.Read[0]["average"])
            except Exception:
                # sra-stat will sometimes put warnings in the XML
                # stream, so we end up with nothing valid to parse.
                # https://github.com/ncbi/sra-tools/issues/192
                logger.error(
                    "Unable to determine index length! Defaulting to small",
                    file=job_context["sra_input_file_path"],
                )
                job_context["index_length_raw"] = -1

    if job_context["index_length_raw"] > 75:
        job_context["index_length"] = "long"
    else:
        job_context["index_length"] = "short"

    # We are trying to determine if the SRA file has single or paired reads. We learned in
    # https://github.com/AlexsLemonade/refinebio/pull/2394 that sometimes sra-stat thinks that a
    # file has paired reads when it really has single reads (like the samples in SRP253951),
    # so we are checking the annotations first then falling back to sra-stat.
    try:
        sample = job_context["sample"]
        for exp in sample.experiments.all():
            for annotation in exp.experimentannotation_set.all():
                if annotation.data.get("library_layout", "").upper() == "PAIRED":
                    job_context["sra_num_reads"] = 2
                    return job_context
                if annotation.data.get("library_layout", "").upper() == "SINGLE":
                    job_context["sra_num_reads"] = 1
                    return job_context

        job_context["sra_num_reads"] = num_reads
    except Exception as e:
        logger.exception(
            "Problem trying to determine library strategy (single/paired)!",
            file=job_context["sra_input_file_path"],
        )
        job_context[
            "job"
        ].failure_reason = "Unable to determine library strategy (single/paired): " + str(e)
        job_context["success"] = False
        return job_context

    if not job_context.get("sra_num_reads", None):
        logger.error(
            "Completely unable to determine library strategy (single/paired)!",
            file=job_context["sra_input_file_path"],
        )
        job_context["job"].failure_reason = "Unable to determine library strategy (single/paired)"
        job_context["success"] = False

    return job_context


def _determine_index_length(job_context: Dict) -> Dict:
    """Determines whether to use the long or short salmon index.

    Adds the key 'index_length' to the job_context with a value of
    'short' if the short index is appropriate or 'long' if the long
    index is appropriate. For more information on index length see the
    _create_index function of the transcriptome_index processor.
    """

    if job_context.get("sra_input_file_path", None):
        return _determine_index_length_sra(job_context)

    logger.debug("Determining index length..")
    total_base_pairs = 0
    number_of_reads = 0
    counter = 1

    if ".gz" == job_context["input_file_path"][-3:]:
        cat = "zcat"
    else:
        cat = "cat"
    # zcat unzips the file provided and dumps the output to STDOUT.
    # It is installed by default in Debian so it should be included
    # in every docker image already.
    with subprocess.Popen(
        [cat, job_context["input_file_path"]], stdout=subprocess.PIPE, universal_newlines=True
    ) as process:
        for line in process.stdout:
            # In the FASTQ file format, there are 4 lines for each
            # read. Three of these contain metadata about the
            # read. The string representing the read itself is found
            # on the second line of each quartet.
            if counter % 4 == 2:
                total_base_pairs += len(line.replace("\n", ""))
                number_of_reads += 1
            counter += 1

    if "input_file_path_2" in job_context:
        if ".gz" == job_context["input_file_path_2"][-3:]:
            cat = "zcat"
        else:
            cat = "cat"
        with subprocess.Popen(
            [cat, job_context["input_file_path_2"]], stdout=subprocess.PIPE, universal_newlines=True
        ) as process:
            for line in process.stdout:
                if counter % 4 == 2:
                    total_base_pairs += len(line.replace("\n", ""))
                    number_of_reads += 1
                counter += 1

    if number_of_reads == 0:
        logger.error(
            "Unable to determine number_of_reads for job.",
            input_file_1=job_context.get("input_file_path"),
            input_file_2=job_context.get("input_file_path_2"),
            job_id=job_context["job"].id,
        )
        job_context["job"].failure_reason = "Unable to determine number_of_reads."
        job_context["job"].no_retry = True
        job_context["success"] = False

        job_context["sample"].is_unable_to_be_processed = True
        job_context["sample"].save()
        return job_context

    index_length_raw = total_base_pairs / number_of_reads

    # Put the raw index length into the job context in a new field for regression testing purposes
    job_context["index_length_raw"] = index_length_raw

    if index_length_raw > 75:
        job_context["index_length"] = "long"
    else:
        job_context["index_length"] = "short"

    return job_context


def _find_or_download_index(job_context: Dict) -> Dict:
    """Finds the appropriate Salmon Index for this experiment.

    Salmon documentation states:

    "If you want to use Salmon in quasi-mapping-based mode, then you
    first have to build an Salmon index for your transcriptome."

    We have used the Data Refinery to build these indices already,
    this function retrieves the location of the correct index for the
    organism and read length and adds it to the job context.
    """
    logger.debug("Fetching and installing index..")

    index_type = "TRANSCRIPTOME_" + job_context["index_length"].upper()
    index_object = (
        OrganismIndex.objects.filter(organism=job_context["organism"], index_type=index_type)
        .order_by("-created_at")
        .first()
    )

    if not index_object:
        logger.error(
            "Could not run Salmon processor without index for organism",
            organism=job_context["organism"],
            processor_job=job_context["job_id"],
            index_type=index_type,
        )
        job_context["job"].no_retry = True
        job_context["job"].failure_reason = "Missing transcriptome index. (" + index_type + ")"
        job_context["success"] = False
        return job_context

    job_context["index_directory"] = index_object.absolute_directory_path

    try:
        # The organism index only needs to be downloaded from S3 once per
        # organism per index length per EBS volume. We don't know if
        # another job has started downloading it yet, started extracting
        # it yet, or already finished and been symlinked to a common
        # location. Therefore check to see if it's happened before we
        # complete each step.
        version_info_path = job_context["index_directory"] + "/versionInfo.json"

        # Something very bad happened and now there are corrupt indexes installed. Nuke 'em.
        if os.path.exists(version_info_path) and (os.path.getsize(version_info_path) == 0):
            logger.error("We have to nuke a zero-valued index directory: " + version_info_path)
            shutil.rmtree(job_context["index_directory"], ignore_errors=True)
            os.makedirs(job_context["index_directory"], exist_ok=True)

        index_tarball = None
        if not os.path.exists(version_info_path):
            # Index is not installed yet, so download it.
            index_file = ComputedFile.objects.filter(result=index_object.result)[0]
            index_tarball = index_file.sync_from_s3(
                path=job_context["work_dir"] + index_file.filename
            )

        index_hard_dir = None
        if not os.path.exists(version_info_path):
            # Index is still not installed yet, so extract it.

            # Create a temporary location to download the index to, which
            # can be symlinked to once extraction is complete.
            index_hard_dir = os.path.join(LOCAL_ROOT_DIR, job_context["job_dir_prefix"]) + "_index/"
            os.makedirs(index_hard_dir)
            with tarfile.open(index_tarball, "r:gz") as index_archive:
                index_archive.extractall(index_hard_dir)

        if not os.path.exists(version_info_path):
            # Index is still not installed yet, so symlink the files we
            # have to where they are expected to reside.
            os.makedirs(job_context["index_directory"], exist_ok=True)

            index_files = [
                "versionInfo.json",
                "duplicate_clusters.tsv",
                "hash.bin",
                "indexing.log",
                "refInfo.json",
                "sa.bin",
                "genes_to_transcripts.txt",
                "header.json",
                "quasi_index.log",
                "rsd.bin",
                "txpInfo.bin",
            ]

            for subfile in index_files:
                # Remove any broken links before symlinking.
                # If a link is broken, os.path.exists returns false because the
                # file pointed to by the link does not exist but
                # os.symlink throws an error because the broken link itself exists.
                # os.path.lexists, however, returns true if a broken link exists
                target_file = job_context["index_directory"] + "/" + subfile
                if os.path.lexists(target_file):
                    os.unlink(target_file)

                os.symlink(index_hard_dir + subfile, target_file)
        elif index_hard_dir:
            # We have failed the race.
            logger.error("We have failed the index extraction race! Removing dead trees.")
            shutil.rmtree(index_hard_dir, ignore_errors=True)
    except Exception as e:
        error_template = "Failed to download or extract transcriptome index for organism {0}: {1}"
        error_message = error_template.format(str(job_context["organism"]), str(e))
        logger.exception(error_message, processor_job=job_context["job_id"])
        job_context["job"].failure_reason = error_message
        job_context["success"] = False

        # Make sure we don't leave an empty index directory lying around.
        shutil.rmtree(index_hard_dir, ignore_errors=True)
        return job_context

    # The index tarball contains a directory named index, so add that
    # to the path where we should put it.
    job_context["genes_to_transcripts_path"] = os.path.join(
        job_context["index_directory"], "genes_to_transcripts.txt"
    )

    job_context["organism_index"] = index_object

    return job_context


def _run_tximport_for_experiment(
    job_context: Dict, experiment: Experiment, quant_files: List[ComputedFile]
) -> Dict:

    # Download all the quant.sf fles for this experiment. Write all
    # their paths to a file so we can pass a path to that to
    # tximport.R rather than having to pass in one argument per
    # sample.
    tximport_path_list_file = job_context["work_dir"] + "tximport_inputs.txt"
    quant_file_paths = {}
    with open(tximport_path_list_file, "w") as input_list:
        for quant_file in quant_files:
            # We create a directory in the work directory for each (quant.sf) file, as
            # tximport assigns column names based on the parent directory name,
            # and we need those names so that we can reassociate withe samples later.
            # ex., a file with absolute_file_path: /processor_job_1/SRR123_output/quant.sf
            # downloads to: /processor_job_2/SRR123_output/quant.sf
            # So the result file has frame "SRR123_output",
            # which we can associate with sample SRR123
            sample_output = (
                job_context["work_dir"] + str(quant_file.absolute_file_path.split("/")[-2]) + "/"
            )
            os.makedirs(sample_output, exist_ok=True)
            quant_work_path = sample_output + quant_file.filename
            quant_file_path = quant_file.get_synced_file_path(path=quant_work_path)
            input_list.write(quant_file_path + "\n")
            quant_file_paths[quant_file_path] = os.stat(quant_file_path).st_size

    rds_filename = "txi_out.RDS"
    rds_file_path = job_context["work_dir"] + rds_filename
    tpm_filename = "gene_lengthScaledTPM.tsv"
    tpm_file_path = job_context["work_dir"] + tpm_filename
    result = ComputationalResult()
    cmd_tokens = [
        "/usr/bin/Rscript",
        "--vanilla",
        "/home/user/data_refinery_workers/processors/tximport.R",
        "--file_list",
        tximport_path_list_file,
        "--gene2txmap",
        job_context["genes_to_transcripts_path"],
        "--rds_file",
        rds_file_path,
        "--tpm_file",
        tpm_file_path,
    ]
    result.time_start = timezone.now()

    logger.debug(
        "Running tximport with: %s",
        str(cmd_tokens),
        processor_job=job_context["job_id"],
        experiment=experiment.id,
    )

    try:
        tximport_result = subprocess.run(cmd_tokens, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        raise utils.ProcessorJobError(
            "Encountered error in R code while running tximport.R: {}".format(str(e)),
            success=False,
            experiment=experiment.id,
        )

    if tximport_result.returncode != 0:
        raise utils.ProcessorJobError(
            "Found non-zero exit code from R code while running tximport.R: {}".format(
                tximport_result.stderr.decode().strip()
            ),
            success=False,
            experiment=experiment.id,
            quant_files=quant_files,
            cmd_tokens=cmd_tokens,
            quant_file_paths=quant_file_paths,
        )

    result.time_end = timezone.now()
    result.commands.append(" ".join(cmd_tokens))
    result.is_ccdl = True
    try:
        processor_key = "TXIMPORT"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        raise utils.ProcessorJobError(
            "Failed to set processor: {}".format(e), success=False, processor_key=processor_key
        )

    result.save()
    job_context["pipeline"].steps.append(result.id)

    rds_file = ComputedFile()
    rds_file.absolute_file_path = rds_file_path
    rds_file.filename = rds_filename
    rds_file.result = result
    rds_file.is_smashable = False
    rds_file.is_qc = False
    rds_file.is_public = True
    rds_file.calculate_sha1()
    rds_file.calculate_size()
    rds_file.save()
    job_context["computed_files"].append(rds_file)

    # Split the tximport result into smashable subfiles
    data = pd.read_csv(tpm_file_path, sep="\t", header=0, index_col=0)
    individual_files = []
    frames = np.split(data, len(data.columns), axis=1)
    for frame in frames:
        # Create sample-specific TPM file.
        sample_file_name = frame.columns.values[0] + "_" + tpm_filename
        frame_path = os.path.join(job_context["work_dir"], sample_file_name)
        frame.to_csv(frame_path, sep="\t", encoding="utf-8")

        # The frame column header is based off of the path, which includes _output.
        sample_accession_code = frame.columns.values[0].replace("_output", "")
        sample = Sample.objects.get(accession_code=sample_accession_code)

        computed_file = ComputedFile()
        computed_file.absolute_file_path = frame_path
        computed_file.filename = sample_file_name
        computed_file.result = result
        computed_file.is_smashable = True
        computed_file.is_qc = False
        computed_file.is_public = True
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.save()
        job_context["computed_files"].append(computed_file)
        job_context["smashable_files"].append(computed_file)

        SampleResultAssociation.objects.get_or_create(sample=sample, result=result)

        # Create association with the RDS file.
        SampleComputedFileAssociation.objects.get_or_create(sample=sample, computed_file=rds_file)

        # Create association with TPM file.
        SampleComputedFileAssociation.objects.get_or_create(
            sample=sample, computed_file=computed_file
        )

        individual_files.append(computed_file)
        job_context["samples"].append(sample)

    # Salmon-processed samples aren't marked as is_processed
    # until they are fully tximported, this value sets that
    # for the end_job function.
    job_context["tximported"] = True
    job_context["individual_files"] = individual_files
    return job_context


def set_tximport_inputs(job_context: Dict) -> Dict:
    """Adds to the job_context a mapping from experiments to a list of their quant files.

    Checks all the experiments which contain a sample from the current
    experiment. If any of them are fully processed (at least with
    salmon-quant) then the return dict will include the experiment
    mapping to a list of paths to the quant.sf file for each sample in
    that experiment.
    """
    experiments = job_context["sample"].experiments.all()

    quantified_experiments = {}
    for experiment in experiments:
        # We only want to consider samples that we actually can run salmon on.
        eligible_samples = experiment.samples.filter(source_database="SRA", technology="RNA-SEQ")
        if not eligible_samples.exists():
            continue

        is_tximport_job = "is_tximport_only" in job_context and job_context["is_tximport_only"]
        salmon_quant_files = get_tximport_inputs_if_eligible(experiment, is_tximport_job)

        if is_tximport_job and salmon_quant_files:
            # If the job is only running tximport, then index_length
            # hasn't been set on the job context because we don't have
            # a raw file to run it on. Therefore pull it from one of
            # the result annotations.

            # Can't just do salmon_quant_results[0] because it's a set.
            index_length = salmon_quant_files[0].result.get_index_length()
            if index_length:
                job_context["index_length"] = index_length
            elif "index_length" not in job_context:
                raise utils.ProcessorJobError(
                    (
                        "Found quant result without an annotation specifying its index length."
                        " Why did this happen?!?"
                    ),
                    success=False,
                    no_retry=True,
                )

        if salmon_quant_files:
            quantified_experiments[experiment] = salmon_quant_files

    job_context["tximport_inputs"] = quantified_experiments

    return job_context


def tximport(job_context: Dict) -> Dict:
    """Run tximport R script based on input quant files and the path
    of genes_to_transcripts.txt.
    """
    tximport_inputs = job_context["tximport_inputs"]

    quantified_experiments = 0
    for experiment, quant_files in tximport_inputs.items():
        job_context = _run_tximport_for_experiment(job_context, experiment, quant_files)
        quantified_experiments += 1

    if (
        quantified_experiments == 0
        and "is_tximport_job" in job_context
        and job_context["is_tximport_only"]
    ):
        raise utils.ProcessorJobError(
            "Tximport job ran on no experiments... Why?!?!?", success=False, no_retry=True
        )

    return job_context


def _run_salmon(job_context: Dict) -> Dict:
    """Runs Salmon Quant."""
    logger.debug("Running Salmon..")

    # Salmon needs to be run differently for different sample types.
    # SRA files also get processed differently as we don't want to use fasterq-dump to extract
    # them to disk.
    if job_context.get("sra_input_file_path", None):

        # Single reads
        if job_context["sra_num_reads"] == 1:

            fifo = "/tmp/barney"
            os.mkfifo(fifo)

            dump_str = "fastq-dump --stdout {input_sra_file} > {fifo} &"
            formatted_dump_command = dump_str.format(
                input_sra_file=job_context["sra_input_file_path"], fifo=fifo
            )
            subprocess.Popen(
                formatted_dump_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )

            command_str = (
                "salmon --no-version-check quant -l A -i {index} "
                "-r {fifo} -p 16 -o {output_directory} --seqBias --dumpEq --writeUnmappedNames"
            )
            formatted_command = command_str.format(
                index=job_context["index_directory"],
                input_sra_file=job_context["sra_input_file_path"],
                fifo=fifo,
                output_directory=job_context["output_directory"],
            )
        # Paired are trickier
        else:

            # Okay, for some reason I can't explain, this only works
            # in the temp directory, otherwise the `tee` part will
            # only output to one or the other of the streams
            # (non-deterministically), but not both. This doesn't
            # appear to happen if the fifos are in tmp.
            alpha = "/tmp/alpha"
            os.mkfifo(alpha)
            beta = "/tmp/beta"
            os.mkfifo(beta)

            dump_str = (
                "fastq-dump --stdout --split-files -I {input_sra_file}"
                "| tee >(grep '@.*\.1\s' -A3 --no-group-separator > {fifo_alpha}) "
                ">(grep '@.*\.2\s' -A3 --no-group-separator > {fifo_beta}) > /dev/null &"
            )
            formatted_dump_command = dump_str.format(
                input_sra_file=job_context["sra_input_file_path"], fifo_alpha=alpha, fifo_beta=beta
            )
            subprocess.Popen(
                formatted_dump_command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            command_str = (
                "salmon --no-version-check quant -l A -i {index} "
                "-1 {fifo_alpha} -2 {fifo_beta} -p 16 -o {output_directory} "
                "--seqBias --dumpEq --writeUnmappedNames"
            )
            formatted_command = command_str.format(
                index=job_context["index_directory"],
                input_sra_file=job_context["sra_input_file_path"],
                fifo_alpha=alpha,
                fifo_beta=beta,
                output_directory=job_context["output_directory"],
            )

    else:
        if "input_file_path_2" in job_context:
            second_read_str = " -2 {}".format(job_context["input_file_path_2"])

            # Rob recommends 16 threads/process, which fits snugly on
            # an x1 at 8GB RAM per Salmon container:

            # (2 threads/core * 16 cores/socket * 64 vCPU) / (1TB/8GB) = ~17
            command_str = (
                "salmon --no-version-check quant -l A --biasSpeedSamp 5 -i {index}"
                " -1 {input_one}{second_read_str} -p 16 -o {output_directory}"
                " --seqBias --gcBias --dumpEq --writeUnmappedNames"
            )

            formatted_command = command_str.format(
                index=job_context["index_directory"],
                input_one=job_context["input_file_path"],
                second_read_str=second_read_str,
                output_directory=job_context["output_directory"],
            )
        else:
            # Related: https://github.com/COMBINE-lab/salmon/issues/83
            command_str = (
                "salmon --no-version-check quant -l A -i {index}"
                " -r {input_one} -p 16 -o {output_directory}"
                " --seqBias --dumpEq --writeUnmappedNames"
            )

            formatted_command = command_str.format(
                index=job_context["index_directory"],
                input_one=job_context["input_file_path"],
                output_directory=job_context["output_directory"],
            )

    logger.debug(
        "Running Salmon Quant using the following shell command: %s",
        formatted_command,
        processor_job=job_context["job_id"],
    )

    # Salmon probably shouldn't take longer than three hours.
    timeout = 60 * 60 * 3
    job_context["time_start"] = timezone.now()
    try:
        completed_command = subprocess.run(
            formatted_command.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        failure_reason = "Salmon timed out because it failed to complete within 3 hours."
        logger.error(
            failure_reason,
            sample_accesion_code=job_context["sample"].accession_code,
            processor_job=job_context["job_id"],
        )
        job_context["job"].failure_reason = failure_reason
        job_context["job"].no_retry = True
        job_context["success"] = False
        return job_context

    job_context["time_end"] = timezone.now()

    if completed_command.returncode == 1:
        stderr = completed_command.stderr.decode().strip()
        error_start = stderr.upper().find("ERROR:")
        error_start = error_start if error_start != -1 else 0
        logger.error(
            "Shell call to salmon failed with error message: %s",
            stderr[error_start:],
            processor_job=job_context["job_id"],
        )

        # If salmon has an error exit code then we don't want to retry it.
        job_context["job"].no_retry = True
        job_context["job"].failure_reason = (
            "Shell call to salmon failed because: " + stderr[error_start:]
        )
        job_context["success"] = False
    else:
        result = ComputationalResult()
        result.commands.append(formatted_command)
        result.time_start = job_context["time_start"]
        result.time_end = job_context["time_end"]
        result.organism_index = job_context["organism_index"]
        result.is_ccdl = True

        try:
            processor_key = "SALMON_QUANT"
            result.processor = utils.find_processor(processor_key)
        except Exception as e:
            return utils.handle_processor_exception(job_context, processor_key, e)

        # Zip up the output of Salmon Quant
        try:
            with tarfile.open(job_context["output_archive"], "w:gz") as tar:
                tar.add(job_context["output_directory"], arcname=os.sep)
        except Exception:
            logger.exception(
                "Exception caught while zipping processed directory %s",
                job_context["output_directory"],
                processor_job=job_context["job_id"],
            )
            failure_template = "Exception caught while zipping processed directory {}"
            job_context["job"].failure_reason = failure_template.format(
                job_context["output_archive"]
            )
            job_context["success"] = False
            return job_context

        salmon_quant_archive = ComputedFile()
        salmon_quant_archive.absolute_file_path = job_context["output_archive"]
        salmon_quant_archive.filename = os.path.split(job_context["output_archive"])[-1]
        salmon_quant_archive.calculate_sha1()
        salmon_quant_archive.calculate_size()
        salmon_quant_archive.is_public = True
        salmon_quant_archive.is_smashable = False
        salmon_quant_archive.is_qc = False

        quant_file = ComputedFile()
        quant_file.s3_bucket = S3_BUCKET_NAME
        timestamp = str(timezone.now().timestamp()).split(".")[0]
        quant_file.s3_key = "quant_files/sample_{0}_{1}_quant.sf".format(
            job_context["sample"].id, timestamp
        )
        quant_file.filename = "quant.sf"
        quant_file.absolute_file_path = job_context["output_directory"] + "quant.sf"
        quant_file.is_public = False
        quant_file.is_smashable = False
        quant_file.is_qc = False
        quant_file.calculate_sha1()
        quant_file.calculate_size()

        # If we're running in the cloud we need to upload the quant.sf
        # file so that it can be used by a job running on any machine
        # to run tximport. We can't use sync_to_s3 though because we
        # have to sync it before we can save the file so it cannot be
        # discovered by other jobs before it is uploaded.
        if settings.RUNNING_IN_CLOUD:
            try:
                S3.upload_file(
                    quant_file.absolute_file_path,
                    quant_file.s3_bucket,
                    quant_file.s3_key,
                    ExtraArgs={"StorageClass": "STANDARD_IA"},
                )
            except Exception as e:
                logger.exception(
                    e, processor_job=job_context["job_id"], sample=job_context["sample"].id
                )
                failure_template = "Exception caught while uploading quantfile to S3: {}"
                job_context["job"].failure_reason = failure_template.format(
                    quant_file.absolute_file_path
                )
                job_context["success"] = False
                return job_context

        # Here select_for_update() is used as a mutex that forces multiple
        # jobs to execute this block of code in serial manner. See:
        # https://docs.djangoproject.com/en/1.11/ref/models/querysets/#select-for-update
        # Theorectically any rows in any table can be locked here, we're
        # locking all existing rows in ComputationalResult table.
        with transaction.atomic():
            ComputationalResult.objects.select_for_update()
            result.save()
            job_context["quant_result"] = result
            quant_file.result = result
            quant_file.save()

            job_context["result"] = result

            job_context["pipeline"].steps.append(result.id)
            SampleResultAssociation.objects.get_or_create(
                sample=job_context["sample"], result=result
            )
            job_context["sample"].most_recent_quant_file = quant_file
            job_context["sample"].save()

            salmon_quant_archive.result = result
            salmon_quant_archive.save()
            job_context["computed_files"].append(salmon_quant_archive)

        kv = ComputationalResultAnnotation()
        kv.data = {
            "index_length": job_context["index_length"],
            "index_length_get": job_context.get("index_length_raw", None),
        }
        kv.result = result
        kv.is_public = True
        kv.save()

        try:
            with open(
                os.path.join(job_context["output_directory"], "lib_format_counts.json")
            ) as lfc_file:
                format_count_data = json.load(lfc_file)
                kv = ComputationalResultAnnotation()
                kv.data = format_count_data
                kv.result = result
                kv.is_public = True
                kv.save()
        except Exception:
            # See: https://github.com/AlexsLemonade/refinebio/issues/1167
            logger.exception(
                "Error parsing Salmon lib_format_counts JSON output!",
                processor_job=job_context["job_id"],
            )

        try:
            with open(
                os.path.join(job_context["output_directory"], "aux_info", "meta_info.json")
            ) as mi_file:
                meta_info = json.load(mi_file)
                kv = ComputationalResultAnnotation()
                kv.data = meta_info
                kv.result = result
                kv.is_public = True
                kv.save()
        except Exception:
            # See: https://github.com/AlexsLemonade/refinebio/issues/1167
            logger.exception(
                "Error parsing Salmon meta_info JSON output!", processor_job=job_context["job_id"]
            )

        job_context["success"] = True

    return job_context


def _run_salmontools(job_context: Dict) -> Dict:
    """Run Salmontools to extract unmapped genes."""

    logger.debug("Running SalmonTools ...")
    unmapped_filename = job_context["output_directory"] + "aux_info/unmapped_names.txt"

    command_str = "salmontools extract-unmapped -u {unmapped_file} -o {output} "
    output_prefix = job_context["salmontools_directory"] + "unmapped_by_salmon"
    command_str = command_str.format(unmapped_file=unmapped_filename, output=output_prefix)
    if "input_file_path_2" in job_context:
        command_str += "-1 {input_1} -2 {input_2}"
        command_str = command_str.format(
            input_1=job_context["input_file_path"], input_2=job_context["input_file_path_2"]
        )
    else:
        command_str += "-r {input_1}"
        command_str = command_str.format(input_1=job_context["input_file_path"])

    start_time = timezone.now()
    logger.debug(
        "Running the following SalmonTools command: %s",
        command_str,
        processor_job=job_context["job_id"],
    )

    completed_command = subprocess.run(
        command_str.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    end_time = timezone.now()

    # As of SalmonTools 0.1.0, completed_command.returncode is always 0,
    # (even if error happens).  completed_command.stderr is not totally
    # reliable either, because it will output the following line even
    # when the execution succeeds:
    #  "There were <N> unmapped reads\n"
    # in which "<N>" is the number of lines in input unmapped_names.txt.
    #
    # As a workaround, we are using a regular expression here to test
    # the status of SalmonTools execution.  Any text in stderr that is
    # not in the above format is treated as error message.
    status_str = completed_command.stderr.decode().strip()
    success_pattern = r"^There were \d+ unmapped reads$"
    if re.match(success_pattern, status_str):
        # Zip up the output of salmontools
        try:
            with tarfile.open(job_context["salmontools_archive"], "w:gz") as tar:
                tar.add(job_context["salmontools_directory"], arcname=os.sep)
        except Exception:
            logger.exception(
                "Exception caught while zipping processed directory %s",
                job_context["salmontools_directory"],
                processor_job=job_context["job_id"],
            )
            failure_template = "Exception caught while zipping salmontools directory {}"
            job_context["job"].failure_reason = failure_template.format(
                job_context["salmontools_archive"]
            )
            job_context["success"] = False
            return job_context

        result = ComputationalResult()
        result.commands.append(command_str)
        result.time_start = start_time
        result.time_end = end_time
        result.is_ccdl = True

        try:
            processor_key = "SALMONTOOLS"
            result.processor = utils.find_processor(processor_key)
        except Exception as e:
            return utils.handle_processor_exception(job_context, processor_key, e)

        result.save()
        job_context["pipeline"].steps.append(result.id)

        assoc = SampleResultAssociation()
        assoc.sample = job_context["sample"]
        assoc.result = result
        assoc.save()

        computed_file = ComputedFile()
        computed_file.filename = job_context["salmontools_archive"].split("/")[-1]
        computed_file.absolute_file_path = job_context["salmontools_archive"]
        computed_file.calculate_sha1()
        computed_file.calculate_size()
        computed_file.is_public = True
        computed_file.is_smashable = False
        computed_file.is_qc = True
        computed_file.result = result
        computed_file.save()
        job_context["computed_files"].append(computed_file)

        assoc = SampleComputedFileAssociation()
        assoc.sample = job_context["sample"]
        assoc.computed_file = computed_file
        assoc.save()

        job_context["result"] = result
        job_context["success"] = True
    else:  # error in salmontools
        logger.error(
            "Shell call to salmontools failed with error message: %s",
            status_str,
            processor_job=job_context["job_id"],
        )
        job_context["job"].failure_reason = (
            "Shell call to salmontools failed because: " + status_str
        )
        job_context["success"] = False

    return job_context


def salmon(job_id: int) -> None:
    """Main processor function for the Salmon Processor.

    Runs salmon quant command line tool, specifying either a long or
    short read length. Also runs Salmontools and Tximport.
    """
    pipeline = Pipeline(name=PipelineEnum.SALMON.value)
    final_context = utils.run_pipeline(
        {"job_id": job_id, "pipeline": pipeline},
        [
            utils.start_job,
            _set_job_prefix,
            _prepare_files,
            _determine_index_length,
            _find_or_download_index,
            _run_salmon,
            set_tximport_inputs,
            tximport,
            _run_salmontools,
            utils.end_job,
        ],
    )
    return final_context

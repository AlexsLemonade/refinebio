from __future__ import absolute_import, unicode_literals
import os
import subprocess
import tarfile
import gzip
import shutil
from typing import Dict
from data_refinery_common.models import File, BatchKeyValue
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


JOB_DIR_PREFIX = "processor_job_"
GENE_TO_TRANSCRIPT_TEMPLATE = "{gene_id}\t{transcript_id}\n"
GENE_TYPE_COLUMN = 2
# Removes each occurrance of ; and "
IDS_CLEANUP_TABLE = str.maketrans({";": None, "\"": None})


def _set_job_prefix(job_context: Dict) -> str:
    job_context["job_dir_prefix"] = JOB_DIR_PREFIX + str(job_context["job_id"])
    return job_context


def _prepare_files(job_context: Dict) -> Dict:
    """Moves the batch's files from the raw directory to the temp directory.

    Also adds the keys "fasta_file_path" and "gtf_file_path" to
    job_context.
    """
    # Transcriptome index processor jobs have only one batch. Each
    # batch has a fasta file and a gtf file
    batch = job_context["batches"][0]
    fasta_file = File.objects.get(batch=batch, raw_format__exact="fa.gz")
    gtf_file = File.objects.get(batch=batch, raw_format__exact="gtf.gz")

    try:
        fasta_file.download_raw_file(job_context["job_dir_prefix"])
        gtf_file.download_raw_file(job_context["job_dir_prefix"])
    except Exception:
        logger.exception("Exception caught while retrieving raw files.",
                         processor_job=job_context["job_id"],
                         batch=batch.id)

        job_context["job"].failure_reason = "Exception caught while retrieving raw files."
        job_context["success"] = False
        return job_context

    # The files are gzipped when download, but rsem-prepare-reference
    # expects them to be unzipped.
    gzipped_fasta_file_path = fasta_file.get_temp_pre_path(job_context["job_dir_prefix"])
    gzipped_gtf_file_path = gtf_file.get_temp_pre_path(job_context["job_dir_prefix"])
    job_context["fasta_file_path"] = gzipped_fasta_file_path.replace(".gz", "")
    job_context["gtf_file_path"] = gzipped_gtf_file_path.replace(".gz", "")

    with gzip.open(gzipped_fasta_file_path, "rb") as gzipped_file, \
            open(job_context["fasta_file_path"], "wb") as gunzipped_file:
        shutil.copyfileobj(gzipped_file, gunzipped_file)

    with gzip.open(gzipped_gtf_file_path, "rb") as gzipped_file, \
            open(job_context["gtf_file_path"], "wb") as gunzipped_file:
        shutil.copyfileobj(gzipped_file, gunzipped_file)

    job_context["fasta_file"] = fasta_file
    job_context["gtf_file"] = gtf_file
    job_context["success"] = True
    return job_context


def _process_gtf(job_context: Dict) -> Dict:
    """Reads in a .gtf file and generates two new files from it.

    The first is a new .gtf file which has all of the pseudogenes
    filtered out of it. The other is a tsv mapping between gene_ids
    and transcript_ids.  Adds the keys "gtf_file_path" and
    "genes_to_transcripts_path" to job_context.
    """
    work_dir = job_context["gtf_file"].get_temp_dir(job_context["job_dir_prefix"])
    filtered_gtf_path = os.path.join(work_dir, "no_pseudogenes.gtf")
    genes_to_transcripts_path = os.path.join(work_dir, "genes_to_transcripts.txt")

    with open(job_context["gtf_file_path"], 'r') as input_gtf, \
            open(filtered_gtf_path, "w") as filtered_gtf, \
            open(genes_to_transcripts_path, "w") as genes_to_transcripts:
        for line in input_gtf:
            # Filter out any lines containing "pseudogene".
            if "pseudogene" in line:
                continue

            filtered_gtf.write(line)
            tab_split_line = line.split("\t")

            # Skip header lines (they're short and contain no tabs):
            if len(tab_split_line) < 2:
                continue

            if tab_split_line[GENE_TYPE_COLUMN] == "transcript":
                ids_column = tab_split_line[-1].translate(IDS_CLEANUP_TABLE)
                split_ids_column = ids_column.split(" ")

                gene_id_index = split_ids_column.index("gene_id") + 1
                gene_id = split_ids_column[gene_id_index]
                transcript_id_index = split_ids_column.index("transcript_id") + 1
                transcript_id = split_ids_column[transcript_id_index]

                genes_to_transcripts.write(
                    GENE_TO_TRANSCRIPT_TEMPLATE.format(gene_id=gene_id,
                                                       transcript_id=transcript_id))

    # Clean up the unfiltered gtf file, which we no longer need.
    os.remove(job_context["gtf_file_path"])
    job_context["gtf_file_path"] = filtered_gtf_path
    job_context["genes_to_transcripts_path"] = genes_to_transcripts_path
    return job_context


def _handle_shell_error(job_context: Dict, stderr: str, command: str) -> None:
    """Logs an error, cleans up, and updates job_context."""
    logger.error("Shell call to {} failed with error message: %s".format(command),
                 stderr,
                 processor_job=job_context["job_id"],
                 batch=job_context["batches"][0])

    job_context["gtf_file"].remove_temp_directory(job_context["job_dir_prefix"])

    # The failure_reason column is only 256 characters wide.
    error_end = 200
    job_context["job"].failure_reason = ("Shell call to {} failed because: ".format(command)
                                         + stderr[:error_end])
    job_context["success"] = False


def _create_index(job_context: Dict) -> Dict:
    """Creates a salmon transcriptome index.

    This index will only be appropriate for use in running salmon on
    transcripts collected from an organism from the same species as
    the Batch. Additionally it will either be a "long" index or a
    "short" index, which means that it will be appropriate for reads
    with a certain range of base pair lengths. The creator of Salmon,
    the esteemed Dr. Rob Patro, has said (via personal communication):

    "For *most* data (i.e. 75bp or longer), the default k should work
    well.  For reads shorter than 75bp ... one should absolutely use a
    shorter k (probably 23 or 21)."
    """
    work_dir = job_context["gtf_file"].get_temp_dir(job_context["job_dir_prefix"])
    job_context["output_dir"] = os.path.join(work_dir, "index")
    rsem_index_dir = os.path.join(work_dir, "rsem_index")
    os.makedirs(rsem_index_dir, exist_ok=True)

    # RSEM takes a prefix path and then all files generated by it will
    # start with that
    rsem_prefix = os.path.join(rsem_index_dir,
                               job_context["fasta_file"].get_base_name())

    rsem_command_string = (
        "rsem-prepare-reference --gtf {gtf_file}"
        " --transcript-to-gene-map {genes_to_transcripts} {fasta_file} {rsem_prefix}"
    )

    rsem_formatted_command = rsem_command_string.format(
        gtf_file=job_context["gtf_file_path"],
        genes_to_transcripts=job_context["genes_to_transcripts_path"],
        fasta_file=job_context["fasta_file_path"],
        rsem_prefix=rsem_prefix
    )

    rsem_completed_command = subprocess.run(rsem_formatted_command.split(),
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)

    if rsem_completed_command.returncode != 0:
        stderr = str(rsem_completed_command.stderr)
        error_start = stderr.find("Error:")
        error_start = error_start if error_start != -1 else 0
        stderr = stderr[error_start:]
        _handle_shell_error(job_context, stderr, "rsem-prepare-reference")
        job_context["success"] = False
        return job_context

    salmon_command_string = ("salmon --no-version-check index -t {rsem_transcripts}"
                             " -i {index_dir} --type quasi -k {kmer_size}")

    kmer_size_property = BatchKeyValue.objects.get(batch_id=job_context["batches"][0].id,
                                                   key="kmer_size")

    # rsem-prepare-reference outputs a transcripts.fa file which needs
    # to be passed into salmon.
    rsem_transcripts = rsem_prefix + ".transcripts.fa"
    salmon_formatted_command = salmon_command_string.format(
        rsem_transcripts=rsem_transcripts,
        index_dir=job_context["output_dir"],
        kmer_size=kmer_size_property.value)

    salmon_completed_command = subprocess.run(salmon_formatted_command.split(),
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)

    if salmon_completed_command.returncode != 0:
        stderr = str(salmon_completed_command.stderr)
        _handle_shell_error(job_context, stderr, "salmon")

    job_context["success"] = True
    return job_context


def _zip_index(job_context: Dict) -> Dict:
    """Zips the index directory into a single .tar.gz file.

    This makes uploading and retrieving the index easier since it will
    only be a single file along with compressing the size of the file
    during storage.
    """
    temp_post_path = job_context["gtf_file"].get_temp_post_path(job_context["job_dir_prefix"])
    try:
        with tarfile.open(temp_post_path, "w:gz") as tar:
            tar.add(job_context["output_dir"],
                    arcname=os.path.basename(job_context["output_dir"]))
    except:
        logger.exception("Exception caught while zipping index directory %s",
                         temp_post_path,
                         processor_job=job_context["job_id"],
                         batch=job_context["batches"][0].id)

        job_context["gtf_file"].remove_temp_directory(job_context["job_dir_prefix"])

        failure_template = "Exception caught while zipping index directory {}"
        job_context["job"].failure_reason = failure_template.format(temp_post_path)
        job_context["success"] = False
        return job_context

    job_context["files_to_upload"] = [job_context["gtf_file"]]
    job_context["success"] = True
    return job_context


def build_transcriptome_index(job_id: int) -> None:
    """The main function for the Transcriptome Index Processor.

    The steps in this process are as follows:
      * First, files are retrieved from Temporary Storage.
      * Next, they are prepared by removing pseudogenes from the gtf file.
      * Next the tool RSEM's prepare-reference is run.
      * Finally the salmon index command is run
    The output of salmon index is a directory which is pushed in full
    to Permanent Storage.
    """
    utils.run_pipeline({"job_id": job_id},
                       [utils.start_job,
                        _set_job_prefix,
                        _prepare_files,
                        _process_gtf,
                        _create_index,
                        _zip_index,
                        utils.upload_processed_files,
                        utils.cleanup_raw_files,
                        utils.end_job])

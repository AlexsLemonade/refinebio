from __future__ import absolute_import, unicode_literals
import os
import subprocess
import tarfile
import gzip
import shutil

from typing import Dict
from django.utils import timezone

from data_refinery_common.models import (
    Organism,
    ProcessorJob,
    OriginalFile,
    ComputationalResult,
    ComputedFile,
    OrganismIndex
)
from data_refinery_workers._version import __version__
from data_refinery_workers.processors import utils
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_env_variable_gracefully

logger = get_and_configure_logger(__name__)


JOB_DIR_PREFIX = "processor_job_"
GENE_TO_TRANSCRIPT_TEMPLATE = "{gene_id}\t{transcript_id}\n"
GENE_TYPE_COLUMN = 2
# Removes each occurrance of ; and "
IDS_CLEANUP_TABLE = str.maketrans({";": None, "\"": None})

ORGANISM_INDEX_BUCKET = get_env_variable_gracefully("S3_TRANSCRIPTOME_INDEX_BUCKET_NAME")


def _compute_paths(job_context: Dict) -> str:
    """Computes the paths for all the directories used/created by this processor.

    Also computes a couple other path-based properties and adds them to the job_context.
    """
    # All files for the job are in the same directory.
    first_file_path = job_context["original_files"][0].absolute_file_path
    job_context["base_file_path"] = '/'.join(first_file_path.split('/')[:-1])
    job_context["work_dir"] = job_context["base_file_path"] + '/' + \
                              JOB_DIR_PREFIX + str(job_context["job_id"])
    os.makedirs(job_context["work_dir"], exist_ok=True)

    job_context["output_dir"] = job_context["work_dir"] + "/" + "index"
    os.makedirs(job_context["output_dir"], exist_ok=True)
    job_context["rsem_index_dir"] = job_context["work_dir"] + "/" + "rsem_index"
    os.makedirs(job_context["rsem_index_dir"], exist_ok=True)


    job_context["organism_name"] = job_context["base_file_path"].split('/')[-1]
    job_context["organism_with_size"] = job_context["organism_name"] + "_" + job_context['length']

    stamp = str(timezone.now().timestamp()).split('.')[0]
    archive_file_name = job_context["organism_name"] + "_" + \
                        job_context['length'] + "_" + stamp + '.tar.gz'
    job_context["computed_archive"] = job_context['base_file_path'] + '/' + archive_file_name

    return job_context


def _prepare_files(job_context: Dict) -> Dict:
    """Adds the keys "fasta_file_path" and "gtf_file_path" to
    job_context.
    """
    original_files = job_context["original_files"]

    for og_file in original_files:
        if "fa.gz" in og_file.source_filename:
            gzipped_fasta_file_path = og_file.absolute_file_path
            job_context["fasta_file"] = og_file
            job_context["fasta_file_path"] = gzipped_fasta_file_path.replace(".gz", "")
            with gzip.open(gzipped_fasta_file_path, "rb") as gzipped_file, \
                    open(job_context["fasta_file_path"], "wb") as gunzipped_file:
                shutil.copyfileobj(gzipped_file, gunzipped_file)
        elif "gtf.gz" in og_file.source_filename:
            gzipped_gtf_file_path = og_file.absolute_file_path
            job_context["gtf_file"] = og_file
            # Gunzip the GTF file into the work directory.
            gunzipped_path = gzipped_gtf_file_path.replace(".gz", "")
            job_context["gtf_file_path"] = os.path.join(job_context["work_dir"],
                                                        gunzipped_path.split("/")[-1])
            with gzip.open(gzipped_gtf_file_path, "rb") as gzipped_file, \
                    open(job_context["gtf_file_path"], "wb") as gunzipped_file:
                shutil.copyfileobj(gzipped_file, gunzipped_file)

    job_context["success"] = True
    return job_context

def _extract_assembly_version(job_context: Dict) -> Dict:
    """Determine the Ensembl assembly version used for this index.

    Ensembl will periodically release updated versions of the
    assemblies which are where the input files for this processor
    comes from.  All divisions other than the main one have identical
    release versions, but we don't know which division these files
    came from so we can't just hit thier API again. Therefore, look at
    the URL we used to get the files because it contains the assembly
    version.

    I'll admit this isn't the most elegant solution, but since the
    transcriptome index's only database model is the OriginalFiles
    until processing is complete, there's no other way to pass this
    information through to this processor without modifying the
    OriginalFile model.

    The URL path we're attempting follows this pattern (defined in the surveyor)
    ftp://ftp.{url_root}/gtf/{species_sub_dir}/{filename_species}.{assembly}.{assembly_version}.gtf.gz
    and we are attempting to extract {assembly_version}.
    """
    original_files = job_context["original_files"]

    for og_file in original_files:
        if ".gtf.gz" in og_file.source_filename:
            extensionless_url = og_file.source_url[:-7]
            version_start_index = extensionless_url.rfind(".") + 1
            job_context["assembly_version"] = extensionless_url[version_start_index:]

    return job_context


def _process_gtf(job_context: Dict) -> Dict:
    """Reads in a .gtf file and generates two new files from it.

    The first is a new .gtf file which has all of the pseudogenes
    filtered out of it. The other is a tsv mapping between gene_ids
    and transcript_ids.  Adds the keys "gtf_file_path" and
    "genes_to_transcripts_path" to job_context.
    """
    filtered_gtf_path = os.path.join(job_context["work_dir"], "no_pseudogenes.gtf")
    # Generate "genes_to_transcripts.txt" in job_context["output_dir"]
    # so that it will be saved later.
    genes_to_transcripts_path = os.path.join(job_context["output_dir"], "genes_to_transcripts.txt")

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
                 processor_job=job_context["job_id"])

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
    # RSEM takes a prefix path and then all files generated by it will
    # start with that

    job_context["salmon_version"] = subprocess.run(['salmon', '--version'],
                                                stderr=subprocess.PIPE,
                                                stdout=subprocess.PIPE).stderr.decode().strip()

    rsem_prefix = os.path.join(job_context["rsem_index_dir"],
                               job_context['base_file_path'].split('/')[-1])

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

    # The extracted fasta file can be large, no need to leave it lying around.
    os.remove(job_context["fasta_file_path"])

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

    # To me, this looks quite ugly! However, I was told we got these values from
    # Rob Patro, the author of Salmon.
    if job_context['length'] is "long":
        job_context['kmer_size'] = "31"
    else:
        job_context['kmer_size'] = "23"

    # rsem-prepare-reference outputs a transcripts.fa file which needs
    # to be passed into salmon.
    rsem_transcripts = rsem_prefix + ".transcripts.fa"
    job_context["salmon_formatted_command"] = salmon_command_string.format(
        rsem_transcripts=rsem_transcripts,
        index_dir=job_context["output_dir"],
        kmer_size=job_context['kmer_size'])

    job_context['time_start'] = timezone.now()
    salmon_completed_command = subprocess.run(job_context["salmon_formatted_command"].split(),
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    job_context['time_end'] = timezone.now()

    if salmon_completed_command.returncode != 0:
        stderr = str(salmon_completed_command.stderr)
        _handle_shell_error(job_context, stderr, "salmon")
        job_context["success"] = False
        return job_context

    job_context["success"] = True
    return job_context


def _zip_index(job_context: Dict) -> Dict:
    """Zips the index directory into a single .tar.gz file.

    This makes uploading and retrieving the index easier since it will
    only be a single file along with compressing the size of the file
    during storage.
    """
    try:
        with tarfile.open(job_context["computed_archive"], "w:gz") as tar:
            tar.add(job_context["output_dir"],
                    arcname=os.path.basename(job_context["output_dir"]))
    except:
        logger.exception("Exception caught while zipping index directory %s",
                         job_context["output_dir"],
                         processor_job=job_context["job_id"]
                         )
        failure_template = "Exception caught while zipping index directory {}"
        job_context["job"].failure_reason = failure_template.format(job_context["output_dir"])
        job_context["success"] = False
        try:
            os.remove(job_context["computed_archive"])
        except Exception:
            pass
        return job_context

    job_context["files_to_upload"] = [job_context["computed_archive"]]
    job_context["success"] = True
    return job_context

def _populate_index_object(job_context: Dict) -> Dict:
    """ """

    result = ComputationalResult()
    result.command_executed = job_context["salmon_formatted_command"]
    result.is_ccdl = True
    result.system_version = __version__
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    result.save()

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context["computed_archive"]
    computed_file.filename = os.path.split(job_context["computed_archive"])[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.result = result
    computed_file.is_smashable = False
    computed_file.is_qc = False
    #computed_file.sync_to_s3(S3_BUCKET_NAME, computed_file.sha1 + "_" + computed_file.filename)
    # TODO here: delete local file after S3 sync
    computed_file.save()

    organism_object = Organism.get_object_for_name(job_context['organism_name'])
    index_object = OrganismIndex()
    index_object.organism = organism_object
    index_object.source_version = job_context["assembly_version"]
    index_object.salmon_version = job_context["salmon_version"]
    index_object.index_type = "TRANSCRIPTOME_" + job_context['length'].upper()
    index_object.result = result
    logger.info("Uploading %s %s to s3" % (job_context['organism_name'], job_context['length']))
    index_object.upload_to_s3(computed_file.absolute_file_path, ORGANISM_INDEX_BUCKET, logger)
    index_object.save()

    job_context['result'] = result
    job_context['computed_file'] = computed_file
    job_context['index'] = index_object

    return job_context

def build_transcriptome_index(job_id: int, length="long") -> None:
    """The main function for the Transcriptome Index Processor.

    The steps in this process are as follows:
      * First, files are retrieved from Temporary Storage.
      * Next, they are prepared by removing pseudogenes from the gtf file.
      * Next the tool RSEM's prepare-reference is run.
      * Finally the salmon index command is run
    The output of salmon index is a directory which is pushed in full
    to Permanent Storage.
    """
    return utils.run_pipeline({"job_id": job_id, "length": length},
                              [utils.start_job,
                               _compute_paths,
                               _prepare_files,
                               _extract_assembly_version,
                               _process_gtf,
                               _create_index,
                               _zip_index,
                               _populate_index_object,
                               #utils.upload_processed_files,
                               #utils.cleanup_raw_files,
                               utils.end_job])

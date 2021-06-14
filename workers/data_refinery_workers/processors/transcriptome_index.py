import gzip
import os
import shutil
import subprocess
import tarfile
from typing import Dict

from django.utils import timezone

from data_refinery_common.job_lookup import PipelineEnum
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Organism,
    OrganismIndex,
    Pipeline,
)
from data_refinery_common.utils import get_env_variable, get_env_variable_gracefully
from data_refinery_workers.processors import utils

logger = get_and_configure_logger(__name__)
JOB_DIR_PREFIX = "processor_job_"
GENE_TO_TRANSCRIPT_TEMPLATE = "{gene_id}\t{transcript_id}\n"
GENE_TYPE_COLUMN = 2
S3_TRANSCRIPTOME_INDEX_BUCKET_NAME = get_env_variable_gracefully(
    "S3_TRANSCRIPTOME_INDEX_BUCKET_NAME", False
)
LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# Removes each occurrance of ; and "
IDS_CLEANUP_TABLE = str.maketrans({";": None, '"': None})


def get_organism_name_from_path(directory_path: str) -> str:
    """Sometimes these files have strain information appended to
    them, so we want to make sure we only use the first two parts of
    a name like Candida_albicans_sc5314_gca_000784635.
    """
    organism_directory = directory_path.split("/")[-1]
    return "_".join(organism_directory.split("_")[:2]).upper()


def _compute_paths(job_context: Dict) -> str:
    """Computes the paths for all the directories used/created by this processor.

    Also computes a couple other path-based properties and adds them to the job_context.
    """
    # All files for the job are in the same directory.
    first_file_path = job_context["original_files"][0].absolute_file_path
    job_context["base_file_path"] = "/".join(first_file_path.split("/")[:-1])
    job_context["work_dir"] = (
        job_context["base_file_path"]
        + "/"
        + job_context["length"].upper()
        + "/"
        + JOB_DIR_PREFIX
        + str(job_context["job_id"])
        + "/"
    )
    try:
        os.makedirs(job_context["work_dir"])
    except Exception as e:
        logger.exception(
            "Could not create work directory for processor job.", job_context=job_context
        )
        job_context["job"].failure_reason = str(e)
        job_context["success"] = False
        return job_context

    job_context["output_dir"] = job_context["work_dir"] + "index/"
    os.makedirs(job_context["output_dir"], exist_ok=True)
    job_context["rsem_index_dir"] = job_context["work_dir"] + "rsem_index/"
    os.makedirs(job_context["rsem_index_dir"], exist_ok=True)

    # We don't actually associated the organism with the job, so we
    # have to determine the organism name from the file
    # path.
    job_context["organism_name"] = get_organism_name_from_path(job_context["base_file_path"])

    stamp = str(timezone.now().timestamp()).split(".")[0]
    archive_file_name = (
        job_context["organism_name"] + "_" + job_context["length"].upper() + "_" + stamp + ".tar.gz"
    )

    job_context["computed_archive"] = job_context["work_dir"] + archive_file_name

    return job_context


def _prepare_files(job_context: Dict) -> Dict:
    """Adds the keys "fasta_file", "fasta_file_path", "gtf_file", and
    "gtf_file_path" to job_context.
    """
    original_files = job_context["original_files"]

    for og_file in original_files:
        if "fa.gz" in og_file.source_filename:
            gzipped_fasta_file_path = og_file.absolute_file_path
            job_context["fasta_file"] = og_file
            new_fasta_filename = gzipped_fasta_file_path.split("/")[-1].replace(".gz", "")
            job_context["fasta_file_path"] = job_context["work_dir"] + new_fasta_filename
            with gzip.open(gzipped_fasta_file_path, "rb") as gzipped_file, open(
                job_context["fasta_file_path"], "wb"
            ) as gunzipped_file:
                shutil.copyfileobj(gzipped_file, gunzipped_file)
        elif "gtf.gz" in og_file.source_filename:
            gzipped_gtf_file_path = og_file.absolute_file_path
            job_context["gtf_file"] = og_file
            # Gunzip the GTF file into the work directory.
            gunzipped_path = gzipped_gtf_file_path.replace(".gz", "")
            job_context["gtf_file_path"] = os.path.join(
                job_context["work_dir"], gunzipped_path.split("/")[-1]
            )
            with gzip.open(gzipped_gtf_file_path, "rb") as gzipped_file, open(
                job_context["gtf_file_path"], "wb"
            ) as gunzipped_file:
                shutil.copyfileobj(gzipped_file, gunzipped_file)

    job_context["success"] = True
    return job_context


def _extract_assembly_information(job_context: Dict) -> Dict:
    """Determine the Ensembl assembly version and name used for this index.

    Ensembl will periodically release updated versions of the
    assemblies which are where the input files for this processor
    comes from.  All divisions other than the main one have identical
    release versions, but we don't know which division these files
    came from so we can't just hit their API again. Therefore, look at
    the URL we used to get the files because it contains the assembly
    version and name.

    I'll admit this isn't the most elegant solution, but since the
    transcriptome index's only database model is the OriginalFiles
    until processing is complete, there's no other way to pass this
    information through to this processor without modifying the
    OriginalFile model.

    The URL path we're attempting follows this pattern (defined in the surveyor)
    ftp://ftp.{url_root}/gtf/{species_sub_dir}/{filename_species}.{assembly_name}.{assembly_version}.gtf.gz
    and we are attempting to extract {assembly_version} and {assembly_name}.
    """
    original_files = job_context["original_files"]

    for og_file in original_files:
        if ".gtf.gz" in og_file.source_filename:
            extensionless_url = og_file.source_url[:-7]
            version_start_index = extensionless_url.rfind(".") + 1
            job_context["assembly_version"] = extensionless_url[version_start_index:]

            # Decrement the index to skip the period.
            versionless_url = extensionless_url[: version_start_index - 1]

            assembly_name_start_index = versionless_url.rfind(".") + 1
            job_context["assembly_name"] = versionless_url[assembly_name_start_index:]

            database_name = "Ensembl"

            # The division name follows right after the first
            # occurence of the assembly version. We don't have a great
            # way to do this since we don't have any special object
            # for organism index until it's created.
            try:
                division_name = versionless_url.split(job_context["assembly_version"])[1][1:].split(
                    "/"
                )[0]

                # If the url is for the main Ensembl then it won't be any of the following
                if division_name not in ["plants", "metazoa", "fungi", "bacteria", "protists"]:
                    database_name = "EnsemblMain"
                else:
                    database_name = "Ensembl" + str.capitalize(division_name)
            except:
                # We use special URLs for tests. I don't love this,
                # but can't think of anything better right now.
                if "data-refinery-test-assets.s3.amazonaws.com" in og_file.source_url:
                    database_name = "EnsemblMain"
                else:
                    job_context[
                        "job"
                    ].failure_reason = "Failed to retrieve/check for division name from url"
                    job_context["success"] = False

            job_context["database_name"] = database_name

    return job_context


def _process_gtf(job_context: Dict) -> Dict:
    """Reads in a .gtf file and generates two new files from it.

    The first is a new .gtf file which has all of the pseudogenes
    filtered out of it. The other is a tsv mapping between gene_ids
    and transcript_ids. Adds the keys "gtf_file_path" and
    "genes_to_transcripts_path" to job_context.
    """
    filtered_gtf_path = os.path.join(job_context["work_dir"], "no_pseudogenes.gtf")
    # Generate "genes_to_transcripts.txt" in job_context["output_dir"]
    # so that it will be included in the computed tarball.
    genes_to_transcripts_path = os.path.join(job_context["output_dir"], "genes_to_transcripts.txt")

    with open(job_context["gtf_file_path"], "r") as input_gtf, open(
        filtered_gtf_path, "w"
    ) as filtered_gtf, open(genes_to_transcripts_path, "w") as genes_to_transcripts:
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
                    GENE_TO_TRANSCRIPT_TEMPLATE.format(gene_id=gene_id, transcript_id=transcript_id)
                )

    # Clean up the unfiltered gtf file, which we no longer need, unless explicitly told otherwise.
    # This setting is used in one of the tests where we download an unzipped gtf
    # from S3 so we don't want to delete it every time.
    if job_context.get("cleanup_gtf", True):
        os.remove(job_context["gtf_file_path"])
    job_context["gtf_file_path"] = filtered_gtf_path
    job_context["genes_to_transcripts_path"] = genes_to_transcripts_path
    return job_context


def _handle_shell_error(job_context: Dict, stderr: str, command: str) -> None:
    """Logs an error, cleans up, and updates job_context."""
    logger.error(
        "Shell call to {} failed with error message: %s".format(command),
        stderr,
        processor_job=job_context["job_id"],
    )

    job_context["job"].failure_reason = "Shell call to {} failed because: ".format(command) + stderr
    job_context["success"] = False


def _create_index(job_context: Dict) -> Dict:
    """Creates a salmon transcriptome index.

    This index will only be appropriate for use in running salmon on
    transcripts collected from an organism from the same species as
    the file. Additionally it will either be a "long" index or a
    "short" index, which means that it will be appropriate for reads
    with a certain range of base pair lengths. The creator of Salmon,
    the esteemed Dr. Rob Patro, has said (via personal communication):

    "For *most* data (i.e. 75bp or longer), the default k should work
    well.  For reads shorter than 75bp ... one should absolutely use a
    shorter k (probably 23 or 21)."
    """
    # RSEM takes a prefix path and then all files generated by it will
    # start with that

    # Version goes to stderr up until the version where it doesn't:
    # https://github.com/COMBINE-lab/salmon/issues/148
    job_context["salmon_version"] = (
        subprocess.run(["salmon", "--version"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        .stdout.decode()
        .strip()
    )

    # TODO: is this providing a filename prefix or a directory or both?
    rsem_prefix = os.path.join(
        job_context["rsem_index_dir"], job_context["base_file_path"].split("/")[-1]
    )

    rsem_command_string = (
        "rsem-prepare-reference --num-threads 16 --gtf {gtf_file}"
        " --transcript-to-gene-map {genes_to_transcripts} {fasta_file} {rsem_prefix}"
    )

    rsem_formatted_command = rsem_command_string.format(
        gtf_file=job_context["gtf_file_path"],
        genes_to_transcripts=job_context["genes_to_transcripts_path"],
        fasta_file=job_context["fasta_file_path"],
        rsem_prefix=rsem_prefix,
    )

    job_context["time_start"] = timezone.now()
    rsem_completed_command = subprocess.run(
        rsem_formatted_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    # The extracted fasta file can be large, no need to leave it lying around.
    os.remove(job_context["fasta_file_path"])

    if rsem_completed_command.returncode != 0:
        stderr = rsem_completed_command.stderr.decode().strip()
        error_start = stderr.upper().find("ERROR:")
        error_start = error_start if error_start != -1 else 0
        stderr = stderr[error_start:]
        _handle_shell_error(job_context, stderr, "rsem-prepare-reference")
        job_context["success"] = False
        return job_context

    salmon_command_string = (
        "salmon --threads=16 --no-version-check index -t {rsem_transcripts}"
        " -i {index_dir} --type quasi -k {kmer_size}"
    )

    # See this function's docstring for more info.
    if job_context["length"] == "long":
        job_context["kmer_size"] = "31"
    else:
        job_context["kmer_size"] = "23"

    # rsem-prepare-reference outputs a transcripts.fa file which needs
    # to be passed into salmon.
    rsem_transcripts = rsem_prefix + ".transcripts.fa"
    job_context["salmon_formatted_command"] = salmon_command_string.format(
        rsem_transcripts=rsem_transcripts,
        index_dir=job_context["output_dir"],
        kmer_size=job_context["kmer_size"],
    )

    salmon_completed_command = subprocess.run(
        job_context["salmon_formatted_command"].split(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    job_context["time_end"] = timezone.now()

    if salmon_completed_command.returncode != 0:
        stderr = salmon_completed_command.stderr.decode().strip()
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
            tar.add(job_context["output_dir"], arcname=os.path.basename(job_context["output_dir"]))
    except Exception:
        logger.exception(
            "Exception caught while zipping index directory %s",
            job_context["output_dir"],
            processor_job=job_context["job_id"],
        )
        failure_template = "Exception caught while zipping index directory {}"
        job_context["job"].failure_reason = failure_template.format(job_context["output_dir"])
        job_context["success"] = False
        try:
            os.remove(job_context["computed_archive"])
        except Exception:
            pass
        return job_context

    job_context["success"] = True
    return job_context


def _populate_index_object(job_context: Dict) -> Dict:
    """ """
    result = ComputationalResult()
    result.commands.append(job_context["salmon_formatted_command"])
    try:
        processor_key = "TX_INDEX"
        result.processor = utils.find_processor(processor_key)
    except Exception as e:
        return utils.handle_processor_exception(job_context, processor_key, e)

    result.is_ccdl = True
    result.time_start = job_context["time_start"]
    result.time_end = job_context["time_end"]
    result.save()
    job_context["pipeline"].steps.append(result.id)

    computed_file = ComputedFile()
    computed_file.absolute_file_path = job_context["computed_archive"]
    computed_file.filename = os.path.split(job_context["computed_archive"])[-1]
    computed_file.calculate_sha1()
    computed_file.calculate_size()
    computed_file.result = result
    computed_file.is_smashable = False
    computed_file.is_qc = False
    computed_file.save()

    organism_object = Organism.get_object_for_name(job_context["organism_name"])
    index_object = OrganismIndex()
    index_object.organism = organism_object
    index_object.database_name = job_context["database_name"]
    index_object.release_version = job_context["assembly_version"]
    index_object.assembly_name = job_context["assembly_name"]
    index_object.salmon_version = job_context["salmon_version"]
    index_object.index_type = "TRANSCRIPTOME_" + job_context["length"].upper()
    # This is where the index will be extracted to.
    index_object.absolute_directory_path = (
        LOCAL_ROOT_DIR
        + "/TRANSCRIPTOME_INDEX/"
        + organism_object.name
        + "/"
        + job_context["length"]
    )
    index_object.result = result

    if S3_TRANSCRIPTOME_INDEX_BUCKET_NAME:
        logger.info(
            "Uploading %s %s to s3",
            job_context["organism_name"],
            job_context["length"],
            processor_job=job_context["job_id"],
        )
        timestamp = str(timezone.now().timestamp()).split(".")[0]
        s3_key = organism_object.name + "_" + index_object.index_type + "_" + timestamp + ".tar.gz"
        sync_result = computed_file.sync_to_s3(S3_TRANSCRIPTOME_INDEX_BUCKET_NAME, s3_key)
        if sync_result:
            computed_file.delete_local_file()
    else:
        logger.warn(
            "S3_TRANSCRIPTOME_INDEX_BUCKET_NAME not configured, therefore %s %s will not be uploaded.",
            job_context["organism_name"],
            job_context["length"],
            processor_job=job_context["job_id"],
        )

    index_object.save()

    # We uploaded the file ourselves since we wanted it to go to a
    # different bucket than end_job would put it in, therefore empty
    # this list so end_job doesn't try to upload it again.
    job_context["computed_files"] = []

    job_context["result"] = result
    job_context["computed_file"] = computed_file
    job_context["index"] = index_object

    # If there's not a long and a short index for this organism yet,
    # don't delete the input.
    # XXX: This will break once we introduce additional versions of these.
    short_indices = OrganismIndex.objects.filter(
        organism=organism_object,
        index_type="TRANSCRIPTOME_SHORT",
        release_version=job_context["assembly_version"],
    )
    long_indices = OrganismIndex.objects.filter(
        organism=organism_object,
        index_type="TRANSCRIPTOME_LONG",
        release_version=job_context["assembly_version"],
    )
    if short_indices.count() < 1 or long_indices.count() < 1:
        # utils.end_job deletes these, so remove them so it doesn't.
        job_context["original_files"] = []

    return job_context


def build_transcriptome_index(job_id: int, length="long", cleanup=None) -> None:
    """The main function for the Transcriptome Index Processor.

    The steps in this process are as follows:
      * First, files are retrieved from Temporary Storage.
      * Next, they are prepared by removing pseudogenes from the gtf file.
      * Next the tool RSEM's prepare-reference is run.
      * Finally the salmon index command is run
    The output of salmon index is a directory which is pushed in full
    to Permanent Storage.
    """
    pipeline = Pipeline(name=PipelineEnum.TX_INDEX.value)
    initial_job_context = {"job_id": job_id, "length": length, "pipeline": pipeline}

    # When running the tests, don't clean up original files so we don't have to
    # keep downloading them.
    if cleanup is not None:
        initial_job_context["cleanup"] = cleanup

    return utils.run_pipeline(
        initial_job_context,
        [
            utils.start_job,
            _compute_paths,
            _prepare_files,
            _extract_assembly_information,
            _process_gtf,
            _create_index,
            _zip_index,
            _populate_index_object,
            utils.end_job,
        ],
    )

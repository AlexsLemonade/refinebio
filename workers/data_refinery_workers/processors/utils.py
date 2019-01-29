import os
import random
import shutil
import signal
import string
import subprocess
import sys
import yaml

from django.conf import settings
from django.utils import timezone
from enum import Enum, unique
from typing import List, Dict, Callable

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    Dataset,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    Pipeline,
    Processor,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    Sample,
)
from data_refinery_common.utils import get_instance_id, get_env_variable, get_env_variable_gracefully


logger = get_and_configure_logger(__name__)
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
DIRNAME = os.path.dirname(os.path.abspath(__file__))
CURRENT_JOB = None


def signal_handler(sig, frame):
    """Signal Handler, works for both SIGTERM and SIGINT"""
    global CURRENT_JOB
    if not CURRENT_JOB:
        sys.exit(0)
    else:
        CURRENT_JOB.start_time = None
        CURRENT_JOB.num_retries = CURRENT_JOB.num_retries - 1
        CURRENT_JOB.save()
        sys.exit(0)


def create_downloader_job(undownloaded_files: OriginalFile) -> bool:
    """Creates a downloader job to download `undownloaded_files`."""
    original_downloader_job = None
    archive_file = None
    for undownloaded_file in undownloaded_files:
        try:
            original_downloader_job = DownloaderJobOriginalFileAssociation.objects.filter(
                original_file=undownloaded_file
            ).latest('id').downloader_job

            # Found the job so we don't need to keep going.
            break
        except DownloaderJobOriginalFileAssociation.DoesNotExist:
            # If there's no association between this file and any
            # downloader jobs, it's most likely because the original
            # file was created after extracting a archive containing
            # multiple files worth of data.
            # The way to handle this is to find that archive and
            # recreate a downloader job FOR THAT. That archive will
            # have the same filename as the file at the end of the
            # 'source_url' field, because that source URL is pointing
            # to the archive we need.
            archive_filename = undownloaded_file.source_url.split("/")[-1]

            # This file or its job might not exist, but we'll wait
            # until we've checked all the files before calling it a
            # failure.
            try:
                archive_file = OriginalFile.objects.filter(filename=archive_filename).first()
                original_downloader_job = DownloaderJobOriginalFileAssociation.objects.filter(
                    original_file=archive_file
                ).latest('id').downloader_job
                # Found the job so we don't need to keep going.
                break
            except:
                pass

    if not original_downloader_job:
        return False

    new_job = DownloaderJob()
    new_job.downloader_task = original_downloader_job.downloader_task
    new_job.accession_code = original_downloader_job.accession_code
    new_job.save()

    if archive_file:
        # If this downloader job is for an archive file, then the
        # files that were passed into this function aren't what need
        # to be directly downloaded, they were extracted out of this
        # archive. The DownloaderJob will re-extract them and set up
        # the associations for the new ProcessorJob.
        # So double check that it still needs downloading because
        # another file that came out of it could have already
        # recreated the DownloaderJob.
        if archive_file.needs_downloading():
            DownloaderJobOriginalFileAssociation.objects.get_or_create(
                downloader_job=new_job,
                original_file=archive_file
            )
    else:
        for original_file in undownloaded_files:
            DownloaderJobOriginalFileAssociation.objects.get_or_create(
                downloader_job=new_job,
                original_file=undownloaded_file
            )

    return True


def prepare_original_files(job_context):
    """ Provision in the Job context for OriginalFile-driven processors
    """
    job = job_context["job"]
    relations = ProcessorJobOriginalFileAssociation.objects.filter(processor_job=job)
    original_files = OriginalFile.objects.filter(id__in=relations.values('original_file_id'))

    if len(original_files) == 0:
        logger.error("No files found.", processor_job=job.id)
        job_context["success"] = False
        job.failure_reason = "No files were found for the job."
        return job_context

    undownloaded_files = set()
    for original_file in original_files:
        if original_file.needs_downloading():
            undownloaded_files.add(original_file)

    if undownloaded_files:
        logger.info(
            ("One or more files found which were missing or not downloaded."
             " Creating downloader jobs for them and deleting this job."),
            processor_job=job.id,
            missing_files=list(undownloaded_files)
        )

        if not create_downloader_job(undownloaded_files):
            failure_reason = "Missing file for processor job but unable to recreate downloader jobs!"
            logger.error(failure_reason, processor_job=job.id)
            job_context["success"] = False
            job.failure_reason = failure_reason
            return job_context

        # If we can't process the data because it's not on the disk we
        # can't mark the job as a success since it obviously didn't
        # succeed. However if we mark it as a failure the job could be
        # retried triggering yet another DownloaderJob to be created
        # to re-download the data. Therefore the best option is to
        # delete this job.
        job.delete()
        job_context["delete_self"] = True

        return job_context

    job_context["original_files"] = original_files
    original_file = job_context['original_files'][0]
    assocs = OriginalFileSampleAssociation.objects.filter(original_file=original_file)
    samples = Sample.objects.filter(id__in=assocs.values('sample_id')).distinct()
    job_context['samples'] = samples
    job_context["computed_files"] = []

    return job_context


def prepare_dataset(job_context):
    """ Provision in the Job context for Dataset-driven processors
    """
    job = job_context["job"]
    relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)

    # This should never be more than one!
    if relations.count() > 1:
        failure_reason = "More than one dataset for processor job!"
        logger.error(failure_reason, processor_job=job.id)
        job_context["success"] = False
        job.failure_reason = failure_reason
        return job_context
    elif relations.count() == 0:
        failure_reason = "No datasets found for processor job!"
        logger.error(failure_reason, processor_job=job.id)
        job_context["success"] = False
        job.failure_reason = failure_reason
        return job_context

    dataset = Dataset.objects.get(id=relations[0].dataset_id)
    dataset.is_processing = True
    dataset.save()

    # Get the samples to smash
    job_context["dataset"] = dataset
    job_context["samples"] = dataset.get_aggregated_samples()
    job_context["experiments"] = dataset.get_experiments()

    # Just in case
    job_context["original_files"] = []
    job_context["computed_files"] = []
    return job_context

def start_job(job_context: Dict):
    """A processor function to start jobs.

    Record in the database that this job is being started and
    retrieves the job's batches from the database and adds them to the
    dictionary passed in with the key 'batches'.
    """
    job = job_context["job"]

    # This job should not have been started.
    if job.start_time is not None and settings.RUNNING_IN_CLOUD:

        if job.success:
            logger.error("ProcessorJob has already completed succesfully - why are we here again? Bad Nomad!",
                job_id=job.id
            )
            job_context["original_files"] = []
            job_context["computed_files"] = []
            job_context['abort'] = True
            return job_context
        if job.success == False:
            logger.error("ProcessorJob has already completed with a fail - why are we here again? Bad Nomad!",
                job_id=job.id
            )
            job_context["original_files"] = []
            job_context["computed_files"] = []
            job_context['abort'] = True
            return job_context

        logger.error("This processor job has already been started!!!", processor_job=job.id)
        raise Exception("processors.start_job called on job %s that has already been started!" % str(job.id))

    # Set up the SIGTERM handler so we can appropriately handle being interrupted.
    # (`docker stop` uses SIGTERM, not SIGINT.)
    # (however, Nomad sends an SIGINT so catch both.)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    job.worker_id = get_instance_id()
    job.worker_version = SYSTEM_VERSION
    job.start_time = timezone.now()
    job.save()

    global CURRENT_JOB
    CURRENT_JOB = job

    logger.debug("Starting processor Job.", processor_job=job.id, pipeline=job.pipeline_applied)

    # Janitors have no requirement
    if job.pipeline_applied not in ["JANITOR"]:
        # Some jobs take OriginalFiles, other take Datasets
        if job.pipeline_applied not in ["SMASHER", "QN_REFERENCE"]:
            job_context = prepare_original_files(job_context)
            if not job_context.get("success", True):
                return job_context
        else:
            job_context = prepare_dataset(job_context)
            if not job_context.get("success", True):
                return job_context
    else:
        # Just in case
        job_context["original_files"] = []
        job_context["computed_files"] = []

    return job_context


def end_job(job_context: Dict, abort=False):
    """A processor function to end jobs.

    Record in the database that this job has completed and that
    the samples have been processed if not aborted.
    """
    job = job_context["job"]

    if "success" in job_context:
        success = job_context["success"]
    else:
        success = True

    if not abort:
        if job_context.get("success", False) and not (job_context["job"].pipeline_applied in ["SMASHER", "QN_REFERENCE"]):

            # Salmon requires the final `tximport` step to be fully `is_processed`.
            mark_as_processed = True
            if (job_context["job"].pipeline_applied == "SALMON" and not job_context.get('tximported', False)):
                mark_as_processed = False

            if mark_as_processed:
                # This handles most of our cases
                unique_experiments = []
                for sample in job_context.get("samples", []):
                    sample.is_processed = True
                    sample.save()
                    if sample.experiments.all().count() > 0:
                        unique_experiments = list(set(unique_experiments + sample.experiments.all()[::1]))

                # Explicitly for the single-salmon scenario
                if 'sample' in job_context:
                    sample = job_context['sample']
                    sample.is_processed = True
                    sample.save()

                for experiment in unique_experiments:
                    experiment.update_num_samples()
                    experiment.save()

    # If we are aborting, it's because we want to do something
    # different, so leave the original files so that "something
    # different" can use them.
    if (success or job.no_retry) and not abort:
        # Cleanup Original Files
        if 'original_files' in job_context:
            for original_file in job_context['original_files']:
                original_file.delete_local_file()

    if success:
        # S3-sync Computed Files
        for computed_file in job_context.get('computed_files', []):
            # Ensure even distribution across S3 servers
            nonce = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(24))
            result = computed_file.sync_to_s3(S3_BUCKET_NAME, nonce + "_" + computed_file.filename)
            if result:
                computed_file.delete_local_file()
    else:
        for computed_file in job_context.get('computed_files', []):
            computed_file.delete_local_file()
            computed_file.delete()

    # If the pipeline includes any steps, save it.
    if 'pipeline' in job_context:
        pipeline = job_context['pipeline']
        if len(pipeline.steps):
            pipeline.save()

    if "work_dir" in job_context and settings.RUNNING_IN_CLOUD:
        shutil.rmtree(job_context["work_dir"], ignore_errors=True)

    job.success = success
    job.end_time = timezone.now()
    job.save()

    if success:
        logger.debug("Processor job completed successfully.",
                    processor_job=job.id,
                    pipeline_applied=job.pipeline_applied)
    else:
        if not job.failure_reason:
            logger.error("Processor job failed without having failure_reason set. FIX ME!!!!!!!!",
                         processor_job=job.id,
                         pipeline_applied=job.pipeline_applied,
                         no_retry=job.no_retry)
        else:
            logger.error("Processor job failed!",
                         processor_job=job.id,
                         pipeline_applied=job.pipeline_applied,
                         no_retry=job.no_retry,
                         failure_reason=job.failure_reason)

    # Return Final Job context so testers can check it
    return job_context


def run_pipeline(start_value: Dict, pipeline: List[Callable]):
    """Runs a pipeline of processor functions.

    start_value must contain a key 'job_id' which is a valid id for a
    ProcessorJob record.

    Each processor fuction must accept a dictionary and return a
    dictionary.

    Any processor function which returns a dictionary containing a key
    of 'success' with a value of False will cause the pipeline to
    terminate with a call to utils.end_job.

    The key 'job' is reserved for the ProcessorJob currently being
    run.  It is required that the dictionary returned by each
    processor function preserve the mapping for 'job' that was passed
    into it.
    """
    job_id = start_value["job_id"]
    try:
        job = ProcessorJob.objects.get(id=job_id)
    except ProcessorJob.DoesNotExist:
        logger.error("Cannot find processor job record.", processor_job=job_id)
        return

    if len(pipeline) == 0:
        logger.error("Empty pipeline specified.",
                     procesor_job=job_id)

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        try:
            last_result = processor(last_result)
        except Exception as e:
            failure_reason = ("Unhandled exception caught while running processor"
                              " function {} in pipeline: ").format(processor.__name__)
            logger.exception(failure_reason,
                             no_retry=job.no_retry,
                             processor_job=job_id)
            last_result["success"] = False
            last_result["job"].failure_reason = failure_reason + str(e)
            return end_job(last_result)

        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor function %s failed. Terminating pipeline.",
                         processor.__name__,
                         processor_job=job_id,
                         failure_reason=last_result["job"].failure_reason)
            return end_job(last_result)

        # We don't want to run end_job at all if the job has deleted
        # itself, which happens if the data for the job was missing.
        if last_result.get("delete_self", False):
            break

        if last_result.get("abort", False):
            return end_job(last_result, abort=True)

    return last_result


@unique
class PipelineEnum(Enum):
    """Hardcoded pipeline names."""

    AGILENT_TWOCOLOR = "Agilent Two Color"
    ARRAY_EXPRESS = "Array Express"
    ILLUMINA = "Illumina"
    NO_OP = "No Op"
    SALMON = "Salmon"
    SMASHER = "Smasher"
    TX_INDEX = "Transcriptome Index"
    QN_REFERENCE = "Quantile Normalization Reference"
    JANITOR = "Janitor"


@unique
class ProcessorEnum(Enum):
    """Hardcoded processor info in each pipeline."""

    # One processor in "Agilent Two Color" pipeline
    AGILENT_TWOCOLOR = {
        "name": "Agilent SCAN TwoColor",
        "docker_img": "not_available_yet",
        "yml_file": "agilent_twocolor.yml"
    }

    # One processor in "Array Express" pipeline
    AFFYMETRIX_SCAN = {
        "name": "Affymetrix SCAN",
        "docker_img": "dr_affymetrix",
        "yml_file": "affymetrix.yml"
    }

    # One processor in "Illumina" pipeline
    ILLUMINA_SCAN = {
        "name": "Illumina SCAN",
        "docker_img": "dr_illumina",
        "yml_file": "illumina.yml"
    }

    # One processor in "No Op" pipeline
    SUBMITTER_PROCESSED = {
        "name": "Submitter-processed",
        "docker_img": "dr_no_op",
        "yml_file": "no_op.yml"
    }

    # Five processors in "Salmon" pipeline
    MULTIQC = {
        "name": "MultiQC",
        "docker_img": "dr_salmon",
        "yml_file": "multiqc.yml"
    }
    FASTERQ_DUMP = {
        "name": "fasterq-dump",
        "docker_img": "dr_salmon",
        "yml_file": "fasterq_dump.yml"
    }
    SALMON_QUANT = {
        "name": "Salmon Quant",
        "docker_img": "dr_salmon",
        "yml_file": "salmon_quant.yml"
    }
    SALMONTOOLS = {
        "name": "Salmontools",
        "docker_img": "dr_salmon",
        "yml_file": "salmontools.yml"
    }
    TXIMPORT = {
        "name": "Tximport",
        "docker_img": "dr_salmon",
        "yml_file": "tximport.yml"
    }

    # One processor in "Smasher" pipeline
    SMASHER = {
        "name": "Smasher",
        "docker_img": "dr_smasher",
        "yml_file": "smasher.yml"
    }

    # One processor in "Transcriptome Index" pipeline
    TX_INDEX = {
        "name": "Transcriptome Index",
        "docker_img": "dr_transcriptome",
        "yml_file": "transcriptome_index.yml"
    }

    QN_REFERENCE = {
        "name": "Quantile Normalization Reference",
        "docker_img": "dr_smasher",
        "yml_file": "qn.yml"
    }

    @classmethod
    def has_key(cls, key):
        """Class method that tells whether a certain key exists."""
        return key in cls.__members__


def get_os_distro():
    """Returns a string of OS distribution.
    Since we are using Docker, this function only considers Linux distribution.
    Alternative files on Linux: /etc/os-release, /etc/lsb-release
    As a matter of fact, "/etc/issue" doesn't exist on Mac OS X.  We can use
    "sw_vers" command to find its OS information.
    A more cross-platform solution is using "platform" module in Python.
    """

    with open('/etc/issue') as distro_fh:
      return distro_fh.readline().strip('\l\n\\n ')


def get_os_pkgs(pkg_list):
    """Returns a dictionay in which each key is the name of an os-lvel
    package and the corresponding value is the package's version.
    This function assumes the package manager is Debian-based (dpkg/apt).
    """

    pkg_info = dict()
    for pkg in pkg_list:
        process_done = subprocess.run(['dpkg-query', '--show', pkg],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        if process_done.returncode:
            raise Exception("OS-level package %s not found: %s" %
                            (pkg, process_done.stderr.decode().strip())
            )

        version = process_done.stdout.decode().strip().split('\t')[-1]
        pkg_info[pkg] = version

    return pkg_info


def get_cmd_lines(cmd_list):
    """Returns a dictionary in which each key is a command string and
    the corresponding value is the command's stripped output.
    """

    cmd_info = dict()
    for cmd in cmd_list:
        process_done = subprocess.run(cmd.split(),
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        if process_done.returncode:
            raise Exception("Failed to run command line '%s': %s" %
                            (cmd, process_done.stderr.decode().strip())
            )

        output_bytes = process_done.stdout
        # Workaround for "salmon --version" and "salmontools --version"
        # commands, whose outputs are both sent to stderr instead of stdout.
        # Alternatively, we could have used "stderr=subprocess.STDOUT" when
        # initializing process_done, but it is probably a better idea to
        # keep stdout and stderr separate.
        base_cmd = cmd.strip().split()[0]
        if base_cmd == 'salmon' or base_cmd == "salmontools":
            output_bytes = process_done.stderr

        cmd_output = output_bytes.decode().strip()
        cmd_info[cmd] = cmd_output

    return cmd_info


def get_pip_pkgs(pkg_list):
    """Returns a dictionary in which each key is the name of a pip-installed
    package and the corresponding value is the package's version.
    Instead of using: `pip show pkg | grep Version | awk '{print $2}'` to get
    each package's version, we save the output of `pip freeze` first, then
    check the version of each input package in pkg_list.  This approach
    launches the subprocess only once and (hopefully) saves some computational
    resource.
    """

    process_done = subprocess.run(['pip', 'freeze'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if process_done.returncode:
        raise Exception("'pip freeze' failed: %s" % process_done.stderr.decode().strip())

    frozen_pkgs = dict()
    for item in process_done.stdout.decode().split():
        name, version = item.split("==")
        frozen_pkgs[name] = version

    pkg_info = dict()
    for pkg in pkg_list:
        try:
            version = frozen_pkgs[pkg]
        except KeyError:
            raise Exception("Pip package not found: %s" % pkg)

        pkg_info[pkg] = version

    return pkg_info


def get_bioc_version():
    """Returns a string that is the version of "Bioconductor" package in R.
    Note that the data frame returned by installed.packages() does NOT include
    a package named "Bioconductor", so we have to launch another R command to
    find "Bioconductor" version.
    """

    r_command = "tools:::.BioC_version_associated_with_R_version()"
    process_done = subprocess.run(['Rscript', '-e', r_command],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if process_done.returncode:
        raise Exception('R command failed to retrieve Bioconductor version: %s' %
                        process_done.stderr.decode().strip()
        )

    version = process_done.stdout.decode().strip().split()[-1]
    version = version[1:-1]  # Remove the leading and trailing non-ascii characters.

    if len(version) == 0:
        raise Exception('Bioconductor not found')

    return version


def get_most_recent_qn_target_for_organism(organism):
    """ Returns a ComputedFile for QN run for an Organism """

    try:
        annotation = ComputationalResultAnnotation.objects.filter(
            data__organism_id=organism.id,
            data__is_qn=True
        ).order_by(
            '-created_at'
        ).first()
        file = annotation.result.computedfile_set.first()
        return file
    except Exception:
        return None


def get_r_pkgs(pkg_list):
    """Returns a dictionary in which each key is the name of a R package
    and the corresponding value is the package's version.
    """

    # Use "Rscript -e <R_commands>" command to get all user-installed R packages.
    r_commands = "packages.df <- as.data.frame(installed.packages()[, c(1, 3:4)]); \
    packages.df <- packages.df[is.na(packages.df$Priority), 1:2, drop=FALSE]; \
    colnames(packages.df) <- NULL; \
    print(packages.df, row.names=FALSE);"

    process_done = subprocess.run(['Rscript', '-e', r_commands],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if process_done.returncode:
        raise Exception('R command failed to retrieves installed packages: %s' %
                        process_done.stderr.decode().strip()
        )

    r_pkgs = dict()
    for item in process_done.stdout.decode().strip().split('\n'):
        name, version = item.strip().split()
        r_pkgs[name] = version

    # "Brainarray" is a collection that consists of 121 ".*ensgprobe" packages.
    # They share the same version number, so we use 'hgu133plus2hsensgprobe'
    # package to report this uniform version.
    ba_proxy_pkg = 'hgu133plus2hsensgprobe'

    pkg_info = dict()
    for pkg in pkg_list:
        if pkg == 'Bioconductor':
            version = get_bioc_version()
        else:
            try:
                version = r_pkgs[pkg] if pkg != "Brainarray" else r_pkgs[ba_proxy_pkg]
            except KeyError:
                raise Exception("R package not found: %s" % pkg)

        pkg_info[pkg] = version

    return pkg_info


def get_checksums(filenames_list):
    """Returns a dictionary in which each key is a file's name and the
    corresponding value is the file's md5 checksum.
    """

    checksums = dict()
    for filename in filenames_list:
        abs_filepath = os.path.join(DIRNAME, filename)
        process_done = subprocess.run(['md5sum', abs_filepath],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        if process_done.returncode:
            raise Exception("md5sum command error:",
                            process_done.stderr.decode().strip())
        checksum_str = process_done.stdout.decode().strip().split()[0]
        checksums[filename] = checksum_str

    return checksums


def get_runtime_env(yml_filename):
    """Reads input YAML filename and returns a dictionary in which each key
    is a category name of runtime environment and the corresponding value
    is an object that includes version information of packages listed in
    that category.
    """

    runtime_env = dict()
    with open(yml_filename) as yml_fh:
        pkgs = yaml.load(yml_fh)
        for pkg_type, pkg_list in pkgs.items():
            if pkg_type == 'os_distribution':
                value = get_os_distro()
            elif pkg_type == 'os_pkg':
                value = get_os_pkgs(pkg_list)
            elif pkg_type == 'cmd_line':
                value = get_cmd_lines(pkg_list)
            elif pkg_type == 'python':
                value = get_pip_pkgs(pkg_list)
            elif pkg_type == 'R':
                value = get_r_pkgs(pkg_list)
            elif pkg_type == 'checksum':
                value = get_checksums(pkg_list)
            else:
                raise Exception("Unknown category in %s: %s" % (yml_filename, pkg_type))

            runtime_env[pkg_type] = value

    return runtime_env


def find_processor(enum_key):
    """Retursn either a newly created Processor record, or the one in
    database that matches the current processor name, version and environment.
    """

    name = ProcessorEnum[enum_key].value['name']
    docker_image = ProcessorEnum[enum_key].value['docker_img']

    # In current implementation, ALWAYS get the runtime environment.
    yml_path = os.path.join(DIRNAME, ProcessorEnum[enum_key].value['yml_file'])
    environment = get_runtime_env(yml_path)
    obj, status = Processor.objects.get_or_create(name=name,
                                                  version=SYSTEM_VERSION,
                                                  docker_image=docker_image,
                                                  environment=environment)
    return obj


def handle_processor_exception(job_context, processor_key, ex):
    err_str = "Failed to set processor: %s" % ex
    logger.error(err_str, job_id=job_context["job"].id, processor=processor_key)
    job_context["job"].failure_reason = err_str
    job_context["success"] = False
    return job_context

import os
import random
import string
import subprocess
import yaml

from enum import Enum, unique
from typing import List, Dict, Callable
from django.utils import timezone
from data_refinery_common.models import (
    ProcessorJob,
    Pipeline,
    Processor,
    Sample,
    OriginalFile,
    Dataset,
    ProcessorJobOriginalFileAssociation,
    ProcessorJobDatasetAssociation,
    OriginalFileSampleAssociation,
    ComputationalResultAnnotation,
    ComputationalResult
)
from data_refinery_workers._version import __version__
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import get_instance_id, get_env_variable

logger = get_and_configure_logger(__name__)
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
DIRNAME = os.path.dirname(os.path.abspath(__file__))


def start_job(job_context: Dict):
    """A processor function to start jobs.

    Record in the database that this job is being started and
    retrieves the job's batches from the database and adds them to the
    dictionary passed in with the key 'batches'.
    """
    job = job_context["job"]

    # This job should not have been started.
    if job.start_time is not None:
        logger.error("This processor job has already been started!!!", processor_job=job.id)
        raise Exception("processors.start_job called on a job that has already been started!")

    job.worker_id = get_instance_id()
    job.worker_version = __version__
    job.start_time = timezone.now()
    job.save()

    logger.info("Starting processor Job.", processor_job=job.id, pipeline=job.pipeline_applied)

    # Some jobs take OriginalFiles, other take Datasets
    if job.pipeline_applied not in ["SMASHER", "QN_REFERENCE"]:
        relations = ProcessorJobOriginalFileAssociation.objects.filter(processor_job=job)
        original_files = OriginalFile.objects.filter(id__in=relations.values('original_file_id'))

        if len(original_files) == 0:
            logger.error("No files found.", processor_job=job.id)
            job_context["success"] = False
            return job_context

        job_context["original_files"] = original_files
        original_file = job_context['original_files'][0]
        assocs = OriginalFileSampleAssociation.objects.filter(original_file=original_file)
        samples = Sample.objects.filter(id__in=assocs.values('sample_id')).distinct()
        job_context['samples'] = samples
        job_context["computed_files"] = []

    else:
        relations = ProcessorJobDatasetAssociation.objects.filter(processor_job=job)

        # This should never be more than one!
        dataset = Dataset.objects.filter(id__in=relations.values('dataset_id')).first()
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
            # This handles most of our cases
            for sample in job_context["samples"]:
                sample.is_processed = True
                sample.save()

            # Explicitly for the single-salmon scenario
            if 'sample' in job_context:
                sample = job_context['sample']
                sample.is_processed = True
                sample.save()

    # S3-sync Original Files
    if 'original_files' in job_context:
        for original_files in job_context['original_files']:
            original_files.delete_local_file()

    # S3-sync Computed Files
    if 'computed_files' in job_context:
        for computed_file in job_context['computed_files']:
            # Ensure even distribution across S3 servers
            nonce = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
            result = computed_file.sync_to_s3(S3_BUCKET_NAME, nonce + "_" + computed_file.filename)
            if result:
                computed_file.delete_local_file()

    # If the pipeline includes any steps, save it.
    if 'pipeline' in job_context:
        pipeline = job_context['pipeline']
        if len(pipeline.steps):
            pipeline.save()

    if "work_dir" in job_context:
        os.rmtree(job_context["work_dir"])

    job.success = success
    job.end_time = timezone.now()
    job.save()

    if success:
        logger.info("Processor job completed successfully.",
                    processor_job=job.id,
                    pipeline_applied=job.pipeline_applied)
    else:
        if not job.failure_reason:
            logger.error("Processor job failed without having failure_reason set. FIX ME!!!!!!!!",
                         processor_job=job.id,
                         pipeline_applied=job.pipeline_applied)
        else:
            logger.info("Processor job failed!",
                        processor_job=job.id,
                        pipeline_applied=job.pipeline_applied,
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
    run.  The key 'batches' is reserved for the Batches that are
    currently being processed.  It is required that the dictionary
    returned by each processor function preserve the mappings for
    'job' and 'batches' that were passed into it.
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
        except Exception:
            logger.exception("Unhandled exception caught while running processor function %s in pipeline",
                             processor.__name__,
                             processor_job=job_id)
            last_result["success"] = False
            return end_job(last_result)

        if "success" in last_result and last_result["success"] is False:
            logger.error("Processor function %s failed. Terminating pipeline.",
                         processor.__name__,
                         processor_job=job_id,
                         failure_reason=last_result["job"].failure_reason)
            return end_job(last_result)

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

    # Four processors in "Salmon" pipeline
    MULTIQC = {
        "name": "MultiQC",
        "docker_img": "dr_salmon",
        "yml_file": "multiqc.yml"
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
        annotation = ComputationalResultAnnotation.objects.filter(data__organism_id=organism.id).order_by('-created_at').first()
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
                                                  version=__version__,
                                                  docker_image=docker_image,
                                                  environment=environment)
    return obj


def handle_processor_exception(job_context, processor_key, ex):
    err_str = "Failed to set processor: %s" % ex
    logger.error(err_str, job_id=job_context["job"].id, processor=processor_key)
    job_context["job"].failure_reason = err_str
    job_context["success"] = False
    return job_context

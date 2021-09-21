import json
import os
import pickle
import random
import shutil
import signal
import string
import subprocess
import sys
from typing import Callable, Dict, List

from django.conf import settings
from django.utils import timezone

import pandas as pd
import yaml

from data_refinery_common.enums import SMASHER_JOB_TYPES, ProcessorEnum, ProcessorPipeline
from data_refinery_common.job_management import create_downloader_job
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import Processor, ProcessorJob, Sample
from data_refinery_common.utils import get_env_variable, get_instance_id

logger = get_and_configure_logger(__name__)
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")
S3_BUCKET_NAME = get_env_variable("S3_BUCKET_NAME", "data-refinery")
S3_QN_TARGET_BUCKET_NAME = get_env_variable("S3_QN_TARGET_BUCKET_NAME", "data-refinery")
DIRNAME = os.path.dirname(os.path.abspath(__file__))
CURRENT_JOB = None


def signal_handler(sig, frame):
    """Signal Handler, works for both SIGTERM and SIGINT"""
    global CURRENT_JOB
    if CURRENT_JOB:
        CURRENT_JOB.success = False
        CURRENT_JOB.end_time = timezone.now()
        CURRENT_JOB.num_retries = CURRENT_JOB.num_retries - 1
        CURRENT_JOB.failure_reason = "Interruped by SIGTERM/SIGINT: " + str(sig)
        CURRENT_JOB.save()

    sys.exit(0)


def prepare_original_files(job_context):
    """Provision in the Job context for OriginalFile-driven processors"""
    job = job_context["job"]
    original_files = job.original_files.all()

    if original_files.count() == 0:
        raise ProcessorJobError("No files were found for the job.", success=False)

    undownloaded_files = set()
    for original_file in original_files:
        if original_file.needs_downloading(job_context["job_id"]):
            if original_file.is_downloaded:
                # If it needs to be downloaded then it's not
                # downloaded and the is_downloaded field should stop
                # lying about that.
                original_file.is_downloaded = False
                original_file.save()

            undownloaded_files.add(original_file)

    if undownloaded_files:
        logger.info(
            (
                "One or more files found which were missing or not downloaded."
                " Creating downloader jobs for them and deleting this job."
            ),
            processor_job=job.id,
            missing_files=list(undownloaded_files),
        )

        was_job_created = create_downloader_job(
            undownloaded_files, processor_job_id=job_context["job_id"], force=True
        )
        if not was_job_created:
            raise ProcessorJobError(
                "Missing file for processor job but unable to recreate downloader jobs!",
                success=False,
            )

        raise ProcessorJobError(
            "We can not process the data because it is not on the disk",
            success=False,
            no_retry=True,  # this job should not be retried again
            abort=True,  # abort the job and don't do anything else
            undownloaded_files=[file.id for file in undownloaded_files],
        )

    job_context["original_files"] = original_files
    first_original_file = original_files.first()
    samples = Sample.objects.filter(original_files=first_original_file)
    job_context["samples"] = samples
    job_context["computed_files"] = []

    return job_context


def prepare_dataset(job_context):
    """Provision in the Job context for Dataset-driven processors"""
    job = job_context["job"]
    job_datasets = job.datasets.all()

    # This should never be more than one!
    if job_datasets.count() > 1:
        raise ProcessorJobError(
            "More than one dataset for processor job!", success=False, no_retry=True
        )
    elif job_datasets.count() == 0:
        raise ProcessorJobError(
            "No datasets found for processor job!", success=False, no_retry=True
        )

    dataset = job_datasets.first()
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

    original_file = job.original_files.first()
    if (
        not job.pipeline_applied == ProcessorPipeline.TXIMPORT.value
        and original_file
        and not original_file.needs_processing(job_context["job_id"])
    ):
        failure_reason = (
            "Sample has a good computed file, it must have been processed, "
            "so it doesn't need to be downloaded! Aborting!"
        )
        logger.error(failure_reason, job_id=job.id, original_file=original_file)
        job_context["original_files"] = []
        job_context["computed_files"] = []
        job_context["abort"] = True
        # Will be saved by end_job.
        job_context["job"].failure_reason = failure_reason
        return job_context

    # Set up the SIGTERM handler so we can appropriately handle being interrupted.
    # (`docker stop` uses SIGTERM, not SIGINT, but better to catch both.)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    # This job should not have been restarted.
    if job.start_time is not None and settings.RUNNING_IN_CLOUD:
        # Let's just log the event and let the job run instead of failing
        # and also reset the endtime and failure reason, since those fields might have been set
        logger.warn(
            "ProcessorJob was restarted by Batch. We do not know why this happened",
            processor_job=job.id,
            success=job.success,
            failure_reason=job.failure_reason,
            start_time=job.start_time,
            end_time=job.end_time,
        )
        job.end_time = None
        job.failure_reason = None

    job.worker_id = get_instance_id()
    job.worker_version = SYSTEM_VERSION
    job.start_time = timezone.now()
    job.save()

    global CURRENT_JOB
    CURRENT_JOB = job

    logger.debug("Starting processor Job.", processor_job=job.id, pipeline=job.pipeline_applied)

    # Janitor jobs don't operate on file objects.
    # Tximport jobs don't need to download the original file, they
    # just need it to know what experiment to process.
    if job.pipeline_applied not in [
        ProcessorPipeline.JANITOR.value,
        ProcessorPipeline.TXIMPORT.value,
    ]:
        # Some jobs take OriginalFiles, other take Datasets
        if ProcessorPipeline[job.pipeline_applied] not in SMASHER_JOB_TYPES:
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
    success = job_context.get("success", True)

    # Upload first so if this fails we can set success = False and let
    # the rest of the function mark it as failed.
    if success:
        # QN reference files go to a special bucket so they can be
        # publicly available.
        if job_context["job"].pipeline_applied == "QN_REFERENCE":
            s3_bucket = S3_QN_TARGET_BUCKET_NAME
        else:
            s3_bucket = S3_BUCKET_NAME

        # S3-sync Computed Files
        for computed_file in job_context.get("computed_files", []):
            # Ensure even distribution across S3 servers
            nonce = "".join(
                random.choice(string.ascii_lowercase + string.digits) for _ in range(24)
            )
            result = computed_file.sync_to_s3(s3_bucket, nonce + "_" + computed_file.filename)

            if result and settings.RUNNING_IN_CLOUD:
                computed_file.delete_local_file()
            elif not result:
                computed_file.delete()
                success = False
                job_context["success"] = False
                job.failure_reason = "Failed to upload computed file."
                break

    if not success:
        for computed_file in job_context.get("computed_files", []):
            computed_file.delete_local_file()
            if computed_file.id:
                computed_file.delete()

        # if the processor job fails mark all datasets as failed
        if ProcessorPipeline[job.pipeline_applied] in SMASHER_JOB_TYPES:
            for dataset in job.datasets.all():
                dataset.failure_reason = job.failure_reason
                dataset.is_processing = False
                dataset.save()

    if not abort:
        if job_context.get("success", False) and not (
            job_context["job"].pipeline_applied
            in [
                ProcessorPipeline.SMASHER.value,
                ProcessorPipeline.QN_REFERENCE.value,
                ProcessorPipeline.CREATE_COMPENDIA.value,
                ProcessorPipeline.CREATE_QUANTPENDIA.value,
                ProcessorPipeline.JANITOR.value,
            ]
        ):
            # Salmon requires the final `tximport` step to be fully `is_processed`.
            mark_as_processed = True
            if job_context["job"].pipeline_applied == "SALMON" and not job_context.get(
                "tximported", False
            ):
                mark_as_processed = False

            if mark_as_processed:
                # This handles most of our cases
                unique_experiments = []
                for sample in job_context.get("samples", []):
                    sample.is_processed = True
                    sample.save()
                    if sample.experiments.all().count() > 0:
                        unique_experiments = list(
                            set(unique_experiments + sample.experiments.all()[::1])
                        )

                # Explicitly for the single-salmon scenario
                if "sample" in job_context:
                    sample = job_context["sample"]
                    sample.is_processed = True
                    sample.save()

                for experiment in unique_experiments:
                    experiment.update_num_samples()

    # If we are aborting, it's because we want to do something
    # different, so leave the original files so that "something
    # different" can use them.
    if (success or job.no_retry) and not abort:
        # Cleanup Original Files unless "cleanup" is set False. This way we
        # don't have to keep downloading files in our tests.
        if "original_files" in job_context and job_context.get("cleanup", True):
            for original_file in job_context["original_files"]:
                if original_file.needs_processing(job.id):
                    original_file.delete_local_file()

    # If the pipeline includes any steps, save it.
    if "pipeline" in job_context:
        pipeline = job_context["pipeline"]
        if len(pipeline.steps):
            pipeline.save()

    if (
        "work_dir" in job_context
        and job_context["job"].pipeline_applied != ProcessorPipeline.CREATE_COMPENDIA.value
        and settings.RUNNING_IN_CLOUD
    ):
        shutil.rmtree(job_context["work_dir"], ignore_errors=True)

    job.abort = abort
    job.success = success
    job.end_time = timezone.now()
    job.save()

    if success:
        logger.debug(
            "Processor job completed successfully.",
            processor_job=job.id,
            pipeline_applied=job.pipeline_applied,
        )
    else:
        if not job.failure_reason:
            logger.error(
                "Processor job failed without having failure_reason set. FIX ME!!!!!!!!",
                processor_job=job.id,
                pipeline_applied=job.pipeline_applied,
                no_retry=job.no_retry,
            )
        else:
            logger.error(
                "Processor job failed!",
                processor_job=job.id,
                pipeline_applied=job.pipeline_applied,
                no_retry=job.no_retry,
                failure_reason=job.failure_reason,
            )

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
        logger.error("Empty pipeline specified.", processor_job=job_id)

    last_result = start_value
    last_result["job"] = job
    for processor in pipeline:
        try:
            last_result = processor(last_result)
        except ProcessorJobError as e:
            e.update_job(job)
            logger.exception(e.failure_reason, processor_job=job.id, **e.context)
            if e.success is False:
                # end_job will use this and set the value
                last_result["success"] = False
            return end_job(last_result, abort=bool(e.abort))
        except Exception as e:
            failure_reason = (
                "Unhandled exception caught while running processor" " function {} in pipeline: "
            ).format(processor.__name__)
            logger.exception(failure_reason, no_retry=job.no_retry, processor_job=job_id)
            last_result["success"] = False
            last_result["job"].failure_reason = failure_reason + str(e)
            return end_job(last_result)

        if "success" in last_result and last_result["success"] is False:
            logger.error(
                "Processor function %s failed. Terminating pipeline.",
                processor.__name__,
                processor_job=job_id,
                failure_reason=last_result["job"].failure_reason,
            )
            return end_job(last_result)

        if last_result.get("abort", False):
            return end_job(last_result, abort=True)

    return last_result


class ProcessorJobError(Exception):
    """General processor job error class."""

    def __init__(
        self, failure_reason, *, success=None, no_retry=None, retried=None, abort=None, **context
    ):
        super(ProcessorJobError, self).__init__(failure_reason)
        self.failure_reason = failure_reason
        self.success = success
        self.no_retry = no_retry
        self.retried = retried
        self.abort = abort
        # additional context to be included when logging
        self.context = context

    def update_job(self, job):
        job.failure_reason = self.failure_reason
        if self.success is not None:
            job.success = self.success
        if self.no_retry is not None:
            job.no_retry = self.no_retry
        if self.retried is not None:
            job.retried = self.retried
        if self.abort is not None:
            job.abort = self.abort
        job.save()

        # also update the failure reason if this is a dataset's processor job
        for dataset in job.datasets.all():
            dataset.failure_reason = self.failure_reason
            dataset.success = False
            dataset.save()


def get_os_distro():
    """Returns a string of OS distribution.
    Since we are using Docker, this function only considers Linux distribution.
    Alternative files on Linux: /etc/os-release, /etc/lsb-release
    As a matter of fact, "/etc/issue" doesn't exist on Mac OS X.  We can use
    "sw_vers" command to find its OS information.
    A more cross-platform solution is using "platform" module in Python.
    """

    with open("/etc/issue") as distro_fh:
        return distro_fh.readline().strip("\l\n\\n ")


def get_os_pkgs(pkg_list):
    """Returns a dictionay in which each key is the name of an os-lvel
    package and the corresponding value is the package's version.
    This function assumes the package manager is Debian-based (dpkg/apt).
    """

    pkg_info = dict()
    for pkg in pkg_list:
        process_done = subprocess.run(
            ["dpkg-query", "--show", pkg], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if process_done.returncode:
            raise Exception(
                "OS-level package %s not found: %s" % (pkg, process_done.stderr.decode().strip())
            )

        version = process_done.stdout.decode().strip().split("\t")[-1]
        pkg_info[pkg] = version

    return pkg_info


def get_cmd_lines(cmd_list):
    """Returns a dictionary in which each key is a command string and
    the corresponding value is the command's stripped output.
    """

    cmd_info = dict()
    for cmd in cmd_list:
        process_done = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process_done.returncode:
            raise Exception(
                "Failed to run command line '%s': %s" % (cmd, process_done.stderr.decode().strip())
            )

        output_bytes = process_done.stdout
        # Workaround for the "salmontools --version"
        # command, whose outputs are sent to stderr instead of stdout.
        # Alternatively, we could have used "stderr=subprocess.STDOUT" when
        # initializing process_done, but it is probably a better idea to
        # keep stdout and stderr separate.
        base_cmd = cmd.strip().split()[0]
        if base_cmd == "salmontools":
            output_bytes = process_done.stderr

        cmd_output = output_bytes.decode().strip()
        cmd_info[cmd] = cmd_output

    return cmd_info


def get_pip_pkgs(pkg_list):
    """Returns a dictionary in which each key is the name of a pip-installed
    package and the corresponding value is the package's version.
    Instead of using: `pip show pkg | grep Version | awk '{print $2}'` to get
    each package's version, we save the output of `pip list --format=json`
    first, then check the version of each input package in the json list.
    This approach launches the subprocess only once and (hopefully) saves
    some computational resource.
    """

    process_done = subprocess.run(
        ["pip", "list", "--format=json"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process_done.returncode:
        raise Exception(
            "'pip list --format=json' failed: %s" % process_done.stderr.decode().strip()
        )

    frozen_pkgs = dict()
    for package in json.loads(process_done.stdout.decode()):
        frozen_pkgs[package["name"]] = package["version"]

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
    process_done = subprocess.run(
        ["Rscript", "-e", r_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process_done.returncode:
        raise Exception(
            "R command failed to retrieve Bioconductor version: %s"
            % process_done.stderr.decode().strip()
        )

    version = process_done.stdout.decode().strip().split()[-1]
    version = version[1:-1]  # Remove the leading and trailing non-ascii characters.

    if len(version) == 0:
        raise Exception("Bioconductor not found")

    return version


def get_r_pkgs(pkg_list):
    """Returns a dictionary in which each key is the name of a R package
    and the corresponding value is the package's version.
    """

    # Use "Rscript -e <R_commands>" command to get all user-installed R packages.
    r_commands = "packages.df <- as.data.frame(installed.packages()[, c(1, 3:4)]); \
    packages.df <- packages.df[is.na(packages.df$Priority), 1:2, drop=FALSE]; \
    colnames(packages.df) <- NULL; \
    print(packages.df, row.names=FALSE);"

    process_done = subprocess.run(
        ["Rscript", "-e", r_commands], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if process_done.returncode:
        raise Exception(
            "R command failed to retrieves installed packages: %s"
            % process_done.stderr.decode().strip()
        )

    r_pkgs = dict()
    for item in process_done.stdout.decode().strip().split("\n"):
        name, version = item.strip().split()
        r_pkgs[name] = version

    # "Brainarray" is a collection that consists of 121 ".*ensgprobe" packages.
    # They share the same version number, so we use 'hgu133plus2hsensgprobe'
    # package to report this uniform version.
    ba_proxy_pkg = "hgu133plus2hsensgprobe"

    pkg_info = dict()
    for pkg in pkg_list:
        if pkg == "Bioconductor":
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
        process_done = subprocess.run(
            ["md5sum", abs_filepath], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if process_done.returncode:
            raise Exception("md5sum command error:", process_done.stderr.decode().strip())
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
        pkgs = yaml.load(yml_fh, Loader=yaml.SafeLoader)
        for pkg_type, pkg_list in pkgs.items():
            if pkg_type == "os_distribution":
                value = get_os_distro()
            elif pkg_type == "os_pkg":
                value = get_os_pkgs(pkg_list)
            elif pkg_type == "cmd_line":
                value = get_cmd_lines(pkg_list)
            elif pkg_type == "python":
                value = get_pip_pkgs(pkg_list)
            elif pkg_type == "R":
                value = get_r_pkgs(pkg_list)
            elif pkg_type == "checksum":
                value = get_checksums(pkg_list)
            else:
                raise Exception("Unknown category in %s: %s" % (yml_filename, pkg_type))

            runtime_env[pkg_type] = value

    return runtime_env


def find_processor(enum_key):
    """Returns either a newly created Processor record, or the one in
    database that matches the current processor name, version and environment.
    """

    name = ProcessorEnum[enum_key].value["name"]
    docker_image = ProcessorEnum[enum_key].value["docker_img"]

    # In current implementation, ALWAYS get the runtime environment.
    yml_path = os.path.join(DIRNAME, ProcessorEnum[enum_key].value["yml_file"])
    environment = get_runtime_env(yml_path)
    obj, status = Processor.objects.get_or_create(
        name=name, version=SYSTEM_VERSION, docker_image=docker_image, environment=environment
    )
    return obj


def handle_processor_exception(job_context, processor_key, ex):
    err_str = "Failed to set processor: %s" % ex
    logger.error(err_str, job_id=job_context["job"].id, processor=processor_key)
    job_context["job"].failure_reason = err_str
    job_context["success"] = False
    return job_context


def cache_keys(*keys, work_dir_key="work_dir"):
    """Decorator to be applied to a pipeline function.
    Returns a new function that calls the original one and caches the given
    keys into the `work_dir`. On the next call it will load those keys (if they
    exist) and add them to the job_context instead of executing the function."""

    def inner(func):
        # generate a unique name for the cache based on the pipeline name
        # and the cached keys
        cache_name = "__".join(list(keys) + [func.__name__])

        def pipeline(job_context):
            cache_path = os.path.join(job_context[work_dir_key], cache_name)
            if os.path.exists(cache_path):
                # cached values exist, load keys from cacke
                try:
                    values = pickle.load(open(cache_path, "rb"))
                    return {**job_context, **values}
                except (OSError, pickle.PickleError):
                    # don't fail if we can't load the cache
                    logger.warning(
                        "Failed to load cached data for pipeline function.",
                        function_name=func.__name__,
                        keys=keys,
                    )

            # execute the actual function
            job_context = func(job_context)

            try:
                # save cached data for the next run
                values = {key: job_context[key] for key in keys}
                pickle.dump(values, open(cache_path, "wb"))
            except Exception:
                # don't fail if we can't save the cache
                logger.warning(
                    "Failed to cache data for pipeline function.",
                    function_name=func.__name__,
                    keys=keys,
                )
            return job_context

        return pipeline

    return inner


def squish_duplicates(data: pd.DataFrame) -> pd.DataFrame:
    """Squish duplicated rows together.
    XXX/TODO: Is mean the appropriate method here?
              We can make this an option in future.
    Discussion here: https://github.com/AlexsLemonade/refinebio/issues/186#issuecomment-395516419
    """
    return data.groupby(data.index, sort=False).mean()

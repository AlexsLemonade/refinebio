import subprocess
import sys
import yaml

from django.core.management.base import BaseCommand
from data_refinery_common.models import Processor
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_workers._version import __version__
from data_refinery_workers.processors.utils import ProcessorEnum

logger = get_and_configure_logger(__name__)


def get_os_distro():
    """Returns a string of OS distribution.
    Since we are using Docker, this function only considers Linux distribution.
    Other alternatives on Linux: /etc/os-release, /etc/lsb-release
    As a matter of fact, "/etc/issue" doesn't exist on Mac OS X.  We can use
    "sw_vers" command to find its OS information.
    A more cross-platform solution is using "platform" module.
    """
    with open('/etc/issue') as distro_fh:
      return distro_fh.readline().strip('\l\n\\n ')


def get_os_pkgs(pkg_list):
    """Returns a dictionary whose key is the name of a os-lvel package
    and value is the package's version. This function assumes the package
    manager is Debian-based (dpkg/apt). It can be a nightmaire to support
    all different package managers on Linux.
    """
    pkg_info = dict()
    for pkg in pkg_list:
        process_done = subprocess.run(['dpkg-query', '--show', pkg],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        if process_done.returncode:
            logger.error("OS-level package %s not found: %s" %
                         (pkg, process_done.stderr.decode().strip())
            )
            sys.exit(1)

        version = process_done.stdout.decode().strip().split('\t')[-1]
        pkg_info[pkg] = version

    return pkg_info


def get_cmd_lines(cmd_list):
    """Returns a dictionary whose key is the name of a command line and
    value is the command's version.  The version is always retrieved by
    "<cmd --version" command.
    """

    cmd_info = dict()
    for cmd in cmd_list:
        args = [cmd, '--version']

        # As of 08/01/2018, "--version" or "-v" option is NOT supported by
        # "rsem-prepare-reference" command, but it is supported by
        # "rsem-calculate-expression", which is another command in RSEM pkg.
        if cmd == "rsem-prepare-reference":
            args[0] = 'rsem-calculate-expression'

        process_done = subprocess.run(args,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        if process_done.returncode:
            logger.error("Command line %s version check failed: %s" %
                         (cmd, process_done.stderr.decode().strip())
            )
            sys.exit(1)

        output_bytes = process_done.stdout
        # This is a workaround for "salmon --version" and "salmontools --version"
        # commands, whose outputs are sent to stderr instead of stdout.
        if not len(output_bytes):
            output_bytes = process_done.stderr

        version = output_bytes.decode().strip().split()[-1]
        cmd_info[cmd] = version

    return cmd_info


def get_pip_pkgs(pkg_list):
    """Returns a dictionary whose key is the name of a pip-installed package
    and value is the package's version.  Instead of using a command like:
      `pip show pkg | grep Version | awk '{print $2}'`
    to get version for each package, we save the output of `pip freeze` as a
    dictionary first, then check the version of packages in pkg_list.
    This approach launches the subprocess only once and (hopefully) saves some
    computational resource.
    """
    process_done = subprocess.run(['pip', 'freeze'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if process_done.returncode:
        logger.error("pip freeze failed: %s" % process_done.stderr.decode().strip())
        sys.exit(1)

    frozen_pkgs = dict()
    for item in process_done.stdout.decode().split():
        name, version = item.split("==")
        frozen_pkgs[name] = version

    pkg_info = dict()
    for pkg in pkg_list:
        try:
            version = frozen_pkgs[pkg]
        except KeyError:
            logger.error("Pip package not found: %s" % pkg)
            sys.exit(1)
        pkg_info[pkg] = version

    return pkg_info


def get_bioc_version():
    """This function returns a string that is the version of "BioConductor"
    package in R.  "BioConductor" package is special in R. It is NOT included
    in the data frame returned by installed.packages().
    """
    r_command = "tools:::.BioC_version_associated_with_R_version()"
    process_done = subprocess.run(['Rscript', '-e', r_command],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    if process_done.returncode:
        logger.error('R command that retrieves BioConductor version failed: %s' %
                     process_done.stderr.decode().strip()
        )
        sys.exit(1)

    version = process_done.stdout.decode().strip().split()[-1]
    version = version[1:-1]  # Remove the leading and trailing non-ascii characters.

    if len(version) == 0:
        logger.error('BioConductor not found')
        sys.exit(1)

    return version


def get_r_pkgs(pkg_list):
    """Returns a dictionary whose key is the name of a R package
    and value is the package's version.
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
        logger.error('R command that retrieves installed packages failed: %s' %
                     process_done.stderr.decode().strip()
        )
        sys.exit(1)

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
        if pkg == 'BioConductor':
            version = get_bioc_version()
        else:
            try:
                version = r_pkgs[pkg] if pkg != "Brainarray" else r_pkgs[ba_proxy_pkg]
            except KeyError:
                logger.error("R package not found: %s" % pkg)
                sys.exit(1)
        pkg_info[pkg] = version

    return pkg_info


def parseYML(yml_filename):
    """Reads input YAML file and returns a dictionary that includes
    package version information.
    """

    version_info = dict()
    with open(yml_filename) as yml_fh:
        pkgs = yaml.load(yml_fh)
        for pkg_type, pkg_list in pkgs.items():
            if pkg_type == 'os_distribution':
                version_info[pkg_type] = get_os_distro()
            elif pkg_type == 'os_pkg':
                version_info[pkg_type] = get_os_pkgs(pkg_list)
            elif pkg_type == 'cmd_line':
                version_info[pkg_type] = get_cmd_lines(pkg_list)
            elif pkg_type == 'python':
                version_info[pkg_type] = get_pip_pkgs(pkg_list)
            elif pkg_type == 'R':
                version_info[pkg_type] = get_r_pkgs(pkg_list)

    return version_info


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--processor",
            type=str,
            help=("The processor's name."))
        parser.add_argument(
            "--yml-file",
            type=str,
            help=("YAML filename that corresponds to the processor."))
        parser.add_argument(
            "--docker-img",
            type=str,
            help=("Name of Docker image repo that hosts the processor."))

    def handle(self, *args, **options):
        if options["processor"] is None:
            logger.error("Input processor name not specified.")
            sys.exit(1)
        if options["yml_file"] is None:
            logger.error("Input YAML filename not specified.")
            sys.exit(1)
        if options["docker_img"] is None:
            logger.error("Input docker image name not specified.")
            sys.exit(1)

        processor_key = options["processor"]
        if not ProcessorEnum.has_key(processor_key):
            logger.error("Invalid processor key: %s" % processor_key)
            sys.exit(1)

        processor_name = ProcessorEnum[processor_key].value
        yml_filename = options["yml_file"]
        try:
            environment = parseYML(yml_filename)
        except Exception as e:
            logger.error("Failed to import YAML file %s: %s" % (yml_filename, e))
            sys.exit(1)

        # If the processor/version pair already exists, update the record;
        # otherwise create a new record.
        Processor.objects.update_or_create(
            name=processor_name,
            version=__version__,
            defaults={'docker_image': options["docker_img"], 'environment': environment}
        )

        sys.exit(0)

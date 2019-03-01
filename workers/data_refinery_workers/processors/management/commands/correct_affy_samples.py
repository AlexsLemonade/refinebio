"""Tries to find and correct mislabelled Affymetrix samples.

See https://github.com/AlexsLemonade/refinebio/issues/968 for more
context, but essentially we had a bug in the surveyor where if a GEO
superseries had an RNA-Seq sample in it, then all the samples would be
labelled as RNA-SEQ, including samples that were actually
Microarray. This command will correct any such Affymetrix samples by
checking their filenames to see if they are .CEL files. If they are
the file is downloaded and we use affy.io to read the header and
determine the platform accesion code and name.
"""

from typing import Dict

from django.core.management.base import BaseCommand
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import *


logger = get_and_configure_logger(__name__)


def _download_file(download_url: str, file_path: str) -> None:
    """ Download a file from GEO.
    """
    # Ensure directory exists
    os.makedirs(file_path.rsplit('/', 1)[0], exist_ok=True)

    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()

        target_file = open(file_path, "wb")
        with closing(urllib.request.urlopen(download_url)) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        urllib.request.urlcleanup()
    except Exception:
        logger.exception("Exception caught while downloading file %s from %s",
                         file_path, download_url)
        raise
    finally:
        target_file.close()

    return True


def _create_ensg_pkg_map() -> Dict:
    """Reads the text file that was generated when installing ensg R
    packages, and returns a map whose keys are chip names and values are
    the corresponding BrainArray ensg package name.
    """
    ensg_pkg_filename = "/home/user/r_ensg_probe_pkgs.txt"
    chip2pkg = dict()
    with open(ensg_pkg_filename) as file_handler:
        for line in file_handler:
            tokens = line.strip("\n").split("\t")
            # tokens[0] is (normalized) chip name,
            # tokens[1] is the package's URL in this format:
            # http://mbni.org/customcdf/<version>/ensg.download/<pkg>_22.0.0.tar.gz
            pkg_name = tokens[1].split("/")[-1].split("_")[0]
            chip2pkg[tokens[0]] = pkg_name

    return chip2pkg


def _determine_brainarray_package(input_file: str) -> Dict:
    """
    """
    try:
        header = ro.r['::']('affyio', 'read.celfile.header')(input_file)
    except RRuntimeError as e:
        error_template = ("Unable to read Affy header in input file {0}"
                          " while running AFFY_TO_PCL due to error: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.info(error_message)
        return None

    # header is a list of vectors. [0][0] contains the package name.
    punctuation_table = str.maketrans(dict.fromkeys(string.punctuation))
    # Normalize header[0][0]
    package_name = header[0][0].translate(punctuation_table).lower()

    # Headers can contain the version "v1" or "v2", which doesn't
    # appear in the brainarray package name. This replacement is
    # brittle, but the list of brainarray packages is relatively short
    # and we can monitor what packages are added to it and modify
    # accordingly. So far "v1" and "v2" are the only known versions
    # which must be accomodated in this way.
    # Related: https://github.com/data-refinery/data-refinery/issues/141
    package_name_without_version = package_name.replace("v1", "").replace("v2", "")
    chip_pkg_map = _create_ensg_pkg_map()
    return chip_pkg_map.get(package_name_without_version, None)


class Command(BaseCommand):
    def handle(self, *args, **options):
        """"""
        for sample in Sample.objects.filter(technology='RNA-SEQ', source_database='GEO'):
            for original_file in sample.original_files.all():
                if original_file.is_CEL_file() :
                    input_file_path = original_file.source_filename
                    download_success = False
                    try:
                        download_success = _download_file(original_file.source_url,
                                                          input_file_path)
                    except:
                        logger.exception("Failed to download %s from %s",
                                    input_file_path,
                                    original_file.source_url)

                    if download_success:
                        brainarray_package = None
                        try:
                            brainarray_package = _determine_brainarray_package(input_file_path)
                        except:
                            logger.info("Failed to detect platform from downloaded file %s.",
                                        input_file_path)

                        if brainarray_package:
                            # If we've detected the platform using affy, then this
                            # is the best source of truth we'll be able to get, so
                            # update the sample to match it.
                            platform_name = get_readable_affymetrix_names()[brainarray_package]

                            sample.platform_accession_code = brainarray_package
                            sample_object.platform_name = platform_name

                    # Regardless of whether we could detect the
                    # platform successfully or not, we definitely know
                    # it's an Affymetrix Microarray because that's the
                    # only one that makes .CEL files.
                    sample.technology = 'MICROARRAY'
                    sample.manufacturer = 'AFFYMETRTIX'
                    sample.save()

                    # If there's other original files associated with
                    # this sample, we don't need them because we
                    # already corrected the platform.
                    break

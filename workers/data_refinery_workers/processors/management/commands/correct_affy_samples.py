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

import os
import shutil
import string
import urllib
from typing import Dict

from django.core.management.base import BaseCommand
from rpy2.rinterface import RRuntimeError
import rpy2.robjects as ro

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import (
    get_env_variable,
    get_readable_affymetrix_names,
)
from data_refinery_common.models import *


logger = get_and_configure_logger(__name__)
CHUNK_SIZE = 1024 * 256 # chunk_size is in bytes


def _download_file(download_url: str, file_path: str) -> None:
    """Download a file from GEO.
    """

    try:
        logger.debug("Downloading file from %s to %s.",
                     download_url,
                     file_path)

        # Ancient unresolved bug. WTF python: https://bugs.python.org/issue27973
        urllib.request.urlcleanup()

        target_file = open(file_path, "wb")
        with urllib.request.urlopen(download_url) as request:
            shutil.copyfileobj(request, target_file, CHUNK_SIZE)

        urllib.request.urlcleanup()
    except Exception:
        logger.exception("Exception caught while downloading file %s from %s",
                         file_path, download_url)
        return False
    finally:
        target_file.close()

    return True


def _determine_brainarray_package(input_file: str) -> Dict:
    """Uses the R package affy.io to read the .CEL file and determine its platform.
    """
    try:
        header = ro.r['::']('affyio', 'read.celfile.header')(input_file)
    except RRuntimeError as e:
        error_template = ("Unable to read Affy header in input file {0}"
                          " while running AFFY_TO_PCL due to error: {1}")
        error_message = error_template.format(input_file, str(e))
        logger.exception(error_message)
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
    return package_name_without_version


class Command(BaseCommand):
    def handle(self, *args, **options):
        """Main function for this command.

        Basically does what is described at the top of this file.
        """
        LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
        work_dir = LOCAL_ROOT_DIR + "/affy_correction/"
        for sample in Sample.objects.filter(technology='RNA-SEQ', source_database='GEO'):
            # Re-create working dir every iteration so we can't run out of disk space.
            shutil.rmtree(work_dir, ignore_errors=True)
            os.makedirs(work_dir, exist_ok=True)
            for original_file in sample.original_files.all():
                if original_file.is_affy_data() :
                    input_file_path = work_dir + original_file.source_filename
                    download_success = _download_file(original_file.source_url,
                                                      input_file_path)

                    if download_success:
                        try:
                            brainarray_package = _determine_brainarray_package(input_file_path)

                            if brainarray_package:
                                logger.info("Determined the package for sample %d is: " + brainarray_package,
                                            sample.id)
                                # If we've detected the platform using affy, then this
                                # is the best source of truth we'll be able to get, so
                                # update the sample to match it.
                                platform_name = get_readable_affymetrix_names()[brainarray_package]

                                sample.platform_accession_code = brainarray_package
                                sample.platform_name = platform_name
                        except:
                            logger.exception("Failed to detect platform from downloaded file %s.",
                                             input_file_path)

                    # Regardless of whether we could detect the
                    # platform successfully or not, we definitely know
                    # it's an Affymetrix Microarray because that's the
                    # only one that makes .CEL files.
                    sample.technology = 'MICROARRAY'
                    sample.manufacturer = 'AFFYMETRIX'
                    sample.save()

                    # If there's other original files associated with
                    # this sample, we don't need them because we
                    # already corrected the platform.
                    break

import re
import string
from typing import Dict

import rpy2.robjects as ro
from rpy2.rinterface import RRuntimeError

from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

ENSG_PKG_FILENAME = "/home/user/r_ensg_probe_pkgs.txt"


def get_platform_from_CEL(cel_file_path: str) -> str:
    """.CEL files have a header which contains platform information.

    This platform information can have some variability to it, but is
    the most reliable way to determine which platform was used to
    generate a sample. We remove this variablility by eliminating
    punctuation and version tags (which aren't part of a platform
    accession).
    """
    try:
        header = ro.r["::"]("affyio", "read.celfile.header")(cel_file_path)
    except RRuntimeError as e:
        error_template = (
            "Unable to read Affy header in input file {0}"
            " while running AFFY_TO_PCL due to error: {1}"
        )
        logger.info(error_template.format(cel_file_path, str(e)))
        raise

    # header is a list of vectors. [0][0] contains the package name.
    # However it contains punctuation which can be variable.
    punctuation_table = str.maketrans(dict.fromkeys(string.punctuation))
    return header[0][0].translate(punctuation_table).lower()

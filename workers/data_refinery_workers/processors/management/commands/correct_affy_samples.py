"""
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import *


# There's no image with both aspera and the R deps needed to detect
# package, so just don't use aspera.
# def _download_file_aspera(download_url: str,
#                           target_file_path: str,
#                           attempt=0) -> bool:
#     """ Download a file to a location using Aspera by shelling out to the `ascp` client. """

#     try:
#         logger.debug("Downloading file from %s to %s via Aspera.",
#                      download_url,
#                      target_file_path)

#         ascp = ".aspera/cli/bin/ascp"
#         key = ".aspera/cli/etc/asperaweb_id_dsa.openssh"
#         url = download_url
#         user = "anonftp"
#         ftp = "ftp-trace.ncbi.nlm.nih.gov"
#         if url.startswith("ftp://"):
#             url = url.replace("ftp://", "")
#         url = url.replace(ftp, "").replace('ftp.ncbi.nlm.nih.gov', '')

#         # Resume level 1, use encryption, unlimited speed
#         command_str = "{} -i {} -k1 -T {}@{}:{} {}".format(ascp, key, user, ftp, url, target_file_path)
#         formatted_command = command_str.format(src=download_url,
#                                                dest=target_file_path)
#         completed_command = subprocess.run(formatted_command.split(),
#                                            stdout=subprocess.PIPE,
#                                            stderr=subprocess.PIPE)

#         # Something went wrong! Else, just fall through to returning True.
#         if completed_command.returncode != 0:

#             stderr = completed_command.stderr.decode().strip()
#             logger.debug("Shell call of `%s` to ascp failed with error message: %s",
#                          formatted_command,
#                          stderr)

#             # Sometimes, GEO fails mysteriously.
#             # Wait a few minutes and try again.
#             if attempt >= 5:
#                 logger.error("All attempts to download accession via ascp failed: %s\nCommand was: %s",
#                              stderr,
#                              formatted_command)
#                 return False
#             else:
#                 time.sleep(30)
#                 return _download_file_aspera(download_url,
#                                              target_file_path,
#                                              attempt + 1
#                                              )
#     except Exception:
#         logger.exception("Exception caught while downloading file from the URL via Aspera: %s",
#                          download_url)
#         return False

#     # If Aspera has given a zero-byte file for some reason, let's back off and retry.
#     if os.path.getsize(target_file_path) < 1:
#         os.remove(target_file_path)
#         if attempt > 5:
#             return False

#         logger.error("Got zero byte ascp download for target, retrying.",
#                      target_url=download_url)
#         time.sleep(10)
#         return _download_file_aspera(download_url,
#                                      target_file_path,
#                                      attempt + 1
#                                      )
#     return True


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
                        logger.info("Failed to download %s from %s",
                                    original_file.source_url,
                                    original_file.source_filename)

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

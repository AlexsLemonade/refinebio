import os
from setuptools import find_packages, setup
import re

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

# code found at: https://stackoverflow.com/a/7071358/6095378
VERSION_FILE = "data_refinery_common/_version.py"
version_line = open(VERSION_FILE, "rt").read()

# The line we're parsing looks something like:
# __version__ = "1.0.1"
# We want to extract the 1.0.1 part with this regex:
VERSION_REGEX = r"^__version__ = ['\"]([^'\"]*)['\"]"
match_object = re.search(VERSION_REGEX, version_line, re.M)
if match_object:
    version_string = match_object.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSION_FILE,))

setup(
    name="data-refinery-common",
    version=version_string,
    packages=find_packages(),
    include_package_data=True,
    install_requires=["django>=1.10.6",
                      "boto3>=1.4.4",
                      "daiquiri>=1.3.0",
                      "billiard>=3.5.0.3"],
    license="BSD License",
    description="Common functionality to be shared between Data Refinery sub-projects.",
    url="https://www.greenelab.com",
    author="Kurt Wheeler",
    author_email="team@greenelab.com",
    classifiers=[
        "Environment :: Web Environment",
        "Framework :: Django",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Ubuntu",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Internet :: WWW/HTTP",
    ],
)

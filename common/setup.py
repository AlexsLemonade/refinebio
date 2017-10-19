import os
from setuptools import find_packages, setup
import re

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

# code found at: https://stackoverflow.com/a/7071358/6095378
VERSIONFILE = "data_refinery_common/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
match_object = re.search(VSRE, verstrline, re.M)
if match_object:
    verstr = match_object.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name="data-refinery-common",
    version=verstr,
    packages=find_packages(),
    include_package_data=True,
    install_requires=["django>=1.10.6",
                      "boto3>=1.4.4",
                      "daiquiri>=1.3.0"],
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

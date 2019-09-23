import os
from setuptools import find_packages, setup
import re

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

VERSION_FILE = "version"
try:
    with open(VERSION_FILE, "rt") as version_file:
        version_string = version_file.read().strip()
except:
    print("Cannot read version to determine System Version."
          " Please create a file common/version containing an up to date System Version.")
    raise

setup(
    name="data-refinery-common",
    version=version_string,
    packages=find_packages(),
    include_package_data=True,
    # These values are based on what is in common/requirements.txt.
    install_requires=["boto3>=1.9.16",
                      "coverage>=4.5.1",
                      "daiquiri>=1.5.0",
                      "django>=2.1.8",
                      "raven>=6.9.0",
                      "requests>=2.10.1",
                      "retrying>=1.3.3",
                      "psycopg2-binary>=2.7.5",
                      "python-nomad>=0.9.0",
                      "rpy2>=3.1.0"
                      ],
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

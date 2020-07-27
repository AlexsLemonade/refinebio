import os

from setuptools import find_packages, setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

VERSION_FILE = "version"
try:
    with open(VERSION_FILE, "rt") as version_file:
        version_string = version_file.read().strip()
except OSError:
    print(
        "Cannot read version to determine System Version."
        " Please create a file workers/version containing an up to date System Version."
    )
    raise


setup(
    name="data-refinery-workers",
    version=version_string,
    packages=find_packages(),
    include_package_data=True,
    license="BSD License",
    description=(
        "Workers to run Downloader and Processor Jobs for" + " the Data Refinery Project."
    ),
    url="https://github.com/data-refinery/data_refinery/tree/master/workers",
    author="Kurt Wheeler",
    author_email="team@greenelab.com",
    classifiers=[
        "Environment :: Web Environment",
        "Framework :: Django",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: Ubuntu",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Topic :: Internet :: WWW/HTTP",
    ],
)

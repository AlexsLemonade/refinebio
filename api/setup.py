import os

from setuptools import find_packages, setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

VERSION_FILE = "version"
try:
    with open(VERSION_FILE, "rt") as version_file:
        version_string = version_file.read().strip()
except:
    print(
        "Cannot read version to determine System Version."
        " Please create a file api/version containing an up to date System Version."
    )
    raise

setup(
    name="data-refinery-api",
    version=version_string,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    license="BSD License",
    description="Refine.bio API",
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

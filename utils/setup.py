import os
from setuptools import find_packages, setup

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

setup(
    name="data-refinery-utils",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["django>=1.10.6",
                      "data_refinery_models>=0.1.0",
                      "boto3>=1.4.4"],
    license="BSD License",
    description="Common util functions for the Data Refinery project",
    url="https://www.greenelab.com",
    author="Kurt Wheeler",
    author_email="kurt.wheeler91@gmail.com",
    classifiers=[
        "Environment :: Web Environment",
        "Framework :: Django",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Internet :: WWW/HTTP",
    ],
)

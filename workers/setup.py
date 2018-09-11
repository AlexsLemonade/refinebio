import os
from setuptools import find_packages, setup
import re


# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

# This is unsafe, but if this isn't set then we should fail!
version_string = os.environ["SYSTEM_VERSION"]

setup(
    name='data-refinery-workers',
    version=version_string,
    packages=find_packages(),
    include_package_data=True,
    license='BSD License',
    description=('Workers to run Downloader and Processor Jobs for'
                 + ' the Data Refinery Project.'),
    url='https://github.com/data-refinery/data_refinery/tree/master/workers',
    author='Kurt Wheeler',
    author_email='team@greenelab.com',
    classifiers=[
        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: Ubuntu',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Internet :: WWW/HTTP',
    ],
)

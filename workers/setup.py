import os
from setuptools import find_packages, setup
import re


# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

# code found at: https://stackoverflow.com/a/7071358/6095378
VERSIONFILE = "data_refinery_workers/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
match_object = re.search(VSRE, verstrline, re.M)
if match_object:
    verstr = match_object.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name='data-refinery-workers',
    version=verstr,
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
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Internet :: WWW/HTTP',
    ],
)

# Refine.bio [![Build Status](https://circleci.com/gh/AlexsLemonade/refinebio/tree/dev.svg?&style=shield)](https://circleci.com/gh/AlexsLemonade/refinebio/) [![codecov](https://codecov.io/gh/AlexsLemonade/refinebio/branch/master/graph/badge.svg)](https://codecov.io/gh/AlexsLemonade/refinebio)

<!-- This section needs to be drastically improved -->
Refine.bio harmonizes petabytes of publicly available biological data into
ready-to-use datasets for cancer researchers and AI/ML scientists.

This README file is about building and running the refine.bio project source code.

If you're interested in simply using the service, you should [go to the website](https://www.refine.bio) or
[read the documentation](http://docs.refine.bio).

Refine.bio currently has four sub-projects contained within this repo:
- [common](./common) Contains code needed by both `foreman` and `workers`.
- [foreman](./foreman) Discovers data to download/process and manages jobs.
- [workers](./workers) Runs Downloader and Processor jobs.
- [infrasctructure](./infrastructure) Manages infrastructure for Refine.bio.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Development](#development)
  - [Git Workflow](#git-workflow)
  - [Installation](#installation)
    - [Automatic](#automatic)
    - [Linux (Manual)](#linux-manual)
    - [Mac (Manual)](#mac-manual)
    - [Virtual Environment](#virtual-environment)
    - [Services](#services)
      - [Postgres](#postgres)
      - [Nomad](#nomad)
    - [Common Dependecies](#common-dependecies)
    - [ElasticSearch](#elasticsearch)
  - [Testing](#testing)
    - [API](#api)
    - [Common](#common)
    - [Foreman](#foreman)
    - [Workers](#workers)
  - [Style](#style)
  - [Gotchas](#gotchas)
    - [R](#r)
- [Running Locally](#running-locally)
  - [API](#api-1)
  - [Surveyor Jobs](#surveyor-jobs)
    - [Sequence Read Archive](#sequence-read-archive)
    - [Ensembl Transcriptome Indices](#ensembl-transcriptome-indices)
  - [Downloader Jobs](#downloader-jobs)
  - [Processor Jobs](#processor-jobs)
  - [Creating Quantile Normalization Reference Targets](#creating-quantile-normalization-reference-targets)
  - [Creating Compendia](#creating-compendia)
  - [Running Tximport Early](#running-tximport-early)
  - [Checking on Local Jobs](#checking-on-local-jobs)
  - [Development Helpers](#development-helpers)
- [Cloud Deployment](#cloud-deployment)
  - [Docker Images](#docker-images)
  - [Git Crypt](#git-crypt)
  - [Terraform](#terraform)
  - [Running Jobs](#running-jobs)
  - [Log Consumption](#log-consumption)
  - [Dumping and Restoring Database Backups](#dumping-and-restoring-database-backups)
  - [Tearing Down](#tearing-down)
- [Support](#support)
- [Meta-README](#meta-readme)
- [License](#license)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Development

### Git Workflow

`refinebio` uses a
[feature branch](http://nvie.com/posts/a-successful-git-branching-model/)
based workflow. New features should be developed on new feature branches, and
pull requests should be sent to the `dev` branch for code review. Merges into
`master` happen at the end of sprints, and tags in `master` correspond to
production releases.

### Installation

To run Refine.bio locally, you will need to have the
prerequisites installed onto your local machine. This will vary depending on
whether you are developing on a Mac or a Linux machine. Linux instructions
have been tested on Ubuntu 16.04 or later, but other Linux distributions
_should_ be able to run the necessary services. Microsoft Windows is currently
unsupported by this project.

__Note__: The install_all.sh script will configure a git pre-commit hook to auto-format your python code.
This will format your code in the same way as the rest of the project, allowing it to pass our linting check.

#### Automatic

The easiest way to run Refine.bio locally is to run `./scripts/install_all.sh`
to install all of the necessary dependencies. As long as you are using a recent
version of Ubuntu or macOS it should work. If you are using another version of
Linux it should still install most of the dependencies as long as you give the
appropriate `INSTALL_CMD` environment variable, but some dependencies may be
named differently in your package manager than in Ubuntu's.

#### Linux (Manual)

The following services will need to be installed:
- Python3 and Pip: `sudo apt-get -y install python3-pip`
- [Docker](https://www.docker.com/community-edition): Be sure to follow the
[post installation steps](https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user)
so Docker does not need sudo permissions.
- [Terraform](https://www.terraform.io/)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries) can be installed on Linux clients with `sudo ./scripts/install_nomad.sh`.
- [pip3](https://pip.pypa.io/en/stable/) can be installed on Linux clients with `sudo apt-get install python3-pip`
- [black](https://black.readthedocs.io/en/stable/) can be installed on Linux clients with `pip3 install black`
- [git-crypt](https://www.agwa.name/projects/git-crypt/)
- [jq](https://stedolan.github.io/jq/)
- [iproute2](https://wiki.linuxfoundation.org/networking/iproute2)
- [shellcheck](https://github.com/koalaman/shellcheck/)

Instructions for installing Docker, Terraform, and Nomad can be found by
following the link for each service. git-crypt, jq, and iproute2 can be installed via
`sudo apt-get install git-crypt jq iproute2 shellcheck`.

#### Mac (Manual)

The following services will need to be installed:
- [Homebrew](https://brew.sh/)
- [Docker for Mac](https://www.docker.com/docker-mac)
- [Terraform](https://www.terraform.io/)
- [Nomad](https://www.nomadproject.io/)
- [git-crypt](https://www.agwa.name/projects/git-crypt/)
- [iproute2mac](https://github.com/brona/iproute2mac)
- [jq](https://stedolan.github.io/jq/)
- [black](https://black.readthedocs.io/en/stable/)
- [shellcheck](https://github.com/koalaman/shellcheck/)

Instructions for installing [Docker](https://www.docker.com/docker-mac) and [Homebrew](https://brew.sh/) can be found by
on their respective homepages.

Once Homebrew is installed, the other required applications can be installed by running: `brew install iproute2mac git-crypt nomad terraform jq black shellcheck`.

Many of the computational processes running are very memory intensive. You will need
to [raise the amount of virtual memory available to
Docker](https://docs.docker.com/docker-for-mac/#advanced) from the default of
2GB to 12GB or 24GB, if possible.

#### Virtual Environment

Run `./scripts/create_virtualenv.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire `refinebio`
repo. Sub-projects each have their own environments managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

#### Services

`refinebio` also depends on Postgres and Nomad. Postgres can be
run in a local Docker container, but Nomad must be run on your
development machine.

##### Postgres

To start a local Postgres server in a Docker container, use:

```bash
./scripts/run_postgres.sh
```

Then, to initialize the database, run:

```bash
./scripts/install_db_docker.sh
```

If you need to access a `psql` shell for inspecting the database, you can use:

```bash
./scripts/run_psql_shell.sh
```

or if you have `psql` installed this command will give you a better shell experience:

```
source scripts/common.sh && PGPASSWORD=mysecretpassword psql -h $(get_docker_db_ip_address) -U postgres -d data_refinery
```

##### Nomad

Similarly, you will need to run a local [Nomad](https://www.nomadproject.io/) service in development mode.

However if you run Linux and you have followed the [installation instructions](#installation), you
can run Nomad with:

```bash
sudo -E ./scripts/run_nomad.sh
```

(_Note:_ This step may take some time because it downloads lots of files.)

Nomad is an orchestration tool which Refine.bio uses to run
`Surveyor`, `Downloader`, `Processor` jobs. Jobs are queued by sending a message to
the Nomad agent, which will then launch a Docker container which runs
the job. If address conflicts emerge, old Docker containers can be purged
with `docker container prune -f`.

#### Common Dependecies

The [common](./common) sub-project contains common code which is
depended upon by the other sub-projects. So before anything else you
should prepare the distribution directory `common/dist` with this
script:

```bash
./scripts/update_models.sh
```

(_Note:_ This step requires the postgres container to be running and initialized.)

Note: there is a small chance this might fail with a `can't stat`, error. If this happens, you have
to manually change permissions on the volumes directory with `sudo chmod -R 740 volumes_postgres`
then re-run the migrations.

#### ElasticSearch

One of the API endpoints is powered by ElasticSearch. ElasticSearch must be running for this functionality to work. A local ElasticSearch instance in a Docker container can be executed with:

```bash
./scripts/run_es.sh
```

And then the ES Indexes (akin to Postgres 'databases') can be created with:

```bash
./scripts/rebuild_es_index.sh
```

### Testing

The end to end tests require a separate Nomad client to be running so
that the tests can be run without interfering with local
development. The second Nomad client can be started with:

```bash
sudo -E ./scripts/run_nomad.sh -e test
```

To run the entire test suite:

```bash
./scripts/run_all_tests.sh
```

(_Note:_ Running all the tests can take some time, especially the first time because it downloads a lot of files.)

You can use the following to get the current status of nomad when running in the test environment.

```
$ source scripts/common.sh
$ set_nomad_test_address
$ nomad status
```

Running the end to end tests is tricky because Nomad's needs to pull images from docker with our code.
We have a docker image registry that runs locally, but you'll need to update it with different images in order to make the code run.
The script `./scripts/prepare_image.sh` can be used to prepare the images before pushing them.

```
$ ./scripts/prepare_image.sh -i downloaders -d localhost:5000
$ docker push localhost:5000/dr_downloaders:latest

$ ./scripts/prepare_image.sh -i no_op -d localhost:5000
$ docker push localhost:5000/dr_no_op:latest
```

That's for the images `downloaders` and `no_op`, the same need to be executed for the other images: `salmon`, `transcriptome`, `illumina` and `affymetrix`.

If you want to debug the status of a specific nomad job you can use:

```
$ nomad status NO_OP_0_2048/dispatch-1567796915-3d7c7c87
$ nomad status f9c1345b
$ nomad logs f9c1345b
```

`f9c1345b` is the allocation id that it's returned in `nomad status`.

These tests will also be run continuously for each commit via CircleCI.

For more granular testing, you can just run the tests for specific parts of the system.

#### API
To just run the API tests:

```bash
./api/run_tests.sh
```

#### Common
To just run the common tests:

```bash
./common/run_tests.sh
```

#### Foreman
To just run the foreman tests:

```bash
./foreman/run_tests.sh
```

#### Workers
To just run the workers tests:

```bash
./workers/run_tests.sh
```

If you only want to run tests with a specific tag, you can do that too. For
example, to run just the salmon tests:

```bash
./workers/run_tests.sh -t salmon
```

All of our worker tests are tagged, generally based on the Docker image required to run them.
Possible values for worker test tags are:
- affymetrix
- agilent
- downloaders
- illumina
- no_op
- qn (short for quantile normalization)
- salmon
- smasher
- transcriptome

### Style

R files in this repo follow
[Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml).
Python Files in this repo follow
[PEP 8](https://www.python.org/dev/peps/pep-0008/). All files (including
Python and R) have a line length limit of 100 characters.

In addition to following pep8, python files must also conform to the formatting style enforced by [black](https://black.readthedocs.io/en/stable/).
`black` is a highly opinionated auto-formatter.
(`black`'s highly opinionated style is a strict sub-set of pep8.)
The easiest way to conform to this style is to run `black . --line-length=100`.
This will auto-format your code.
Running the `./scripts/install_all.sh` script will install a pre-commit git hook that will run this formatter on every commit you make locally. Under the hood this uses [pre-commit](https://pre-commit.com/), which you can also install directly by running `pip3 install pre-commit & pre-commit install`. Then, if you want to run `pre-commit` without making a git commit, you can use `pre-commit run --all-files`.
To install `black` see the [installation instructions](#installation).
Any Pull Requests that do not conform to the style enforced by `black` will be flagged by our continous integration and will not be accepted until that check passes.

All user-facing scripts have been linted with `shellcheck` for common
warnings and POSIX-correctness. If a script is user-facing, it should ideally
be POSIX-compliant and have the extension `.sh`, but if bashisms are necessary
it should have the extension `.bash`. To install `shellcheck`, you can run
`apt-get install shellcheck` or `brew install shellcheck`. Then, you can lint
scripts with `shellcheck FILE`.



### Gotchas

During development, you make encounter some occasional strangeness. Here's
some things to watch out for:

  - Since we use multiple Docker instances, don't forget to `./scripts/update_models`
  - If builds are failing, increase the size of Docker's memory allocation. (Mac only.)
  - If Docker images are failing mysteriously during creation, it may
be the result of Docker's `Docker.qcow2` or `Docker.raw` file filling. You
can prune old images with `docker system prune -a`.
  - If it's killed abruptly, the containerized Postgres images can be
  left in an unrecoverable state. Annoying.

#### R

We have created some utilities to help us keep R stable, reliable, and from periodically causing build errors related to version incompatibilites.
The primary goal of these is to pin the version for every R package that we have.
The R package `devtools` is useful for this, but in order to be able to install a specific version of it, we've created the R script `common/install_devtools.R`.

There is another gotcha to be aware of should you ever need to modify versions of R or its packages.
In Dockerfiles for images that need the R language, we install apt packages that look like `r-base-core=3.4.2-1xenial1`.
It's unclear why the version for these is so weird, but it was determined by visiting the package list here: https://cran.revolutionanalytics.com/bin/linux/ubuntu/xenial/
If it needs to be updated then a version should be selected from that list.

Additionally there are two apt packages, r-base and r-base-core, which seem to be very similar except that r-base-core is slimmed down some by not including some additional packages.
For a while we were using r-base, but we switched to r-base-core when we pinned the version of the R language because the r-base package caused an apt error.

## Running Locally

Once you've built the `common/dist` directory and have
the Nomad and Postgres services running, you're ready to run the API and/or jobs.
There are three kinds of jobs within Refine.bio.

### API

The API can be run with:

```bash
./api/serve.sh
```

### Surveyor Jobs

Surveyor Jobs discover samples to download/process along with
recording metadata about the samples. A Surveyor Job should queue
`Downloader Jobs` to download the data it discovers.

The Surveyor can be run with the `./foreman/run_surveyor.sh`
script. The first argument to this script is the type of Surveyor Job
to run, which will always be `survey_all`.

Details on these expected arguments can be viewed by
running:

```bash
./foreman/run_surveyor.sh survey_all -h
```

The Surveyor can accept a single accession code from any of the source
data repositories (e.g., Sequencing Read Archive,
 ArrayExpress, Gene Expression Omnibus):

```bash
./foreman/run_surveyor.sh survey_all --accession <ACCESSION_CODE>
```

Example for a GEO experiment:

```bash
./foreman/run_surveyor.sh survey_all --accession GSE85217
```

Example for an ArrayExpress experiment:

```bash
./foreman/run_surveyor.sh survey_all --accession E-MTAB-3050 # AFFY
./foreman/run_surveyor.sh survey_all --accession E-GEOD-3303 # NO_OP
```

Transcriptome indices are a bit special.
For species within the "main" Ensembl division, the species name can be provided like so:

```bash
./foreman/run_surveyor.sh survey_all --accession "Homo sapiens"
```

However for species that are in other divisions, the division must follow the species name after a comma like so:

```bash
./foreman/run_surveyor.sh survey_all --accession "Caenorhabditis elegans, EnsemblMetazoa"
```
The possible divisions that can be specified are:
* Ensembl (this is the "main" division and is the default)
* EnsemblPlants
* EnsemblFungi
* EnsemblBacteria
* EnsemblProtists
* EnsemblMetazoa

If you are unsure what division a species falls into, unfortunately the only way to tell is go to check ensembl.com.
(Although googling the species name + "ensembl" may work pretty well.)

You can also supply a newline-deliminated file to `survey_all` which will
dispatch survey jobs based on accession codes like so:

```bash
./foreman/run_surveyor.sh survey_all --file MY_BIG_LIST_OF_CODES.txt
```

The main foreman job loop can be started with:

```bash
./foreman/run_surveyor.sh retry_jobs
```

This must actually be running for jobs to move forward through the pipeline.

#### Sequence Read Archive

When surveying SRA, you can supply _either_ run accession codes (e.g.,
codes beginning in `SRR`, `DRR`, or `ERR`) or study accession codes
(`SRP`, `DRP`, `ERP`).

Run example (single read):

```bash
./foreman/run_surveyor.sh survey_all --accession DRR002116
```

Run example (paired read):

```bash
./foreman/run_surveyor.sh survey_all --accession SRR6718414
```

Study example:

```bash
./foreman/run_surveyor.sh survey_all --accession ERP006872
```

#### Ensembl Transcriptome Indices

Building transcriptome indices used for quantifying RNA-seq data requires
us to retrieve genome information from
[Ensembl](http://ensemblgenomes.org/). The Surveyor expects a species'
scientific name in the main Ensembl division as the accession:

```bash
./foreman/run_surveyor.sh survey_all --accession "Homo Sapiens"
```

TODO: Update once this supports organisms from multiple Ensembl divisions

### Downloader Jobs

Downloader Jobs will be queued automatically when `Surveyor Jobs`
discover new samples. However, if you just want to queue a `Downloader Job`
yourself rather than having the Surveyor do it for you, you can use the `./workers/tester.sh`
script:
```bash
./workers/tester.sh run_downloader_job --job-name=<EXTERNAL_SOURCE> --job-id=<JOB_ID>
```

For example:
```bash
./workers/tester.sh run_downloader_job --job-name=SRA --job-id=12345
```

or

```bash
./workers/tester.sh run_downloader_job --job-name=ARRAY_EXPRESS --job-id=1
```

Or for more information run:
```bash
./workers/tester.sh -h
```

### Processor Jobs

Processor Jobs will be queued automatically by successful `Downloader Jobs`.
However, if you just want to run a `Processor Job` without yourself without having
a `Downloader Job` do it for you, the following command will do so:

```bash
./workers/tester.sh -i <IMAGE_NAME> run_processor_job --job-name=<JOB_NAME> --job-id=<JOB_ID>
```

For example
```bash
./workers/tester.sh -i affymetrix run_processor_job --job-name=AFFY_TO_PCL --job-id=54321
```

or

```bash
./workers/tester.sh -i no_op run_processor_job --job-name=NO_OP --job-id=1
```

or

```bash
./workers/tester.sh -i salmon run_processor_job --job-name=SALMON --job-id=1
```

or

```bash
./workers/tester.sh -i transcriptome run_processor_job --job-name=TRANSCRIPTOME_INDEX_LONG --job-id=1
```

Or for more information run:
```bash
./workers/tester.sh -h
```

### Creating Quantile Normalization Reference Targets

If you want to quantile normalize combined outputs, you'll first need to create a reference target for a given organism. This can be done in a production environment with the following:

```bash
nomad job dispatch -meta ORGANISM=DANIO_RERIO CREATE_QN_TARGET
```

To create QN targets for all organisms, do so with the dispatcher:

```bash
nomad job dispatch QN_DISPATCHER
```

This will at some point move to the foreman and then it will take a list of organisms to create QN targets for.

### Creating Compendia

Creating species-wide compendia for a given species can be done in a production environment by running the following on the Foreman instance:

```bash
./run_management_command.sh create_compendia --organisms=DANIO_RERIO --svd-algorithm=ARPACK
```

or for a list of organisms:

```bash
./run_management_command.sh create_compendia --organisms=DANIO_RERIO,HOMO_SAPIENS --svd-algorithm=ARPACK
```

or for all organisms with sufficient data:

```bash
./run_management_command.sh create_compendia --svd-algorithm=ARPACK
```

Alternatively a compendium can be created which only includes quant.sf files by using the create_quantpentida command:

```
./run_management_command.sh create_quantpendia --organisms=DANIO_RERIO
```

Compendia jobs run on the smasher instance.
However they require a very large amount of RAM to be able to complete.
Our smasher instance does not generally have enough RAM to be able to run them, so if you need to run a smasher job you should temporarily increase the size of the smasher instance.
This can be done by changing the terraform variable `smasher_instance_type` which can be found in `infrastructure/variables.tf`.
Select an AWS instance type that has enough RAM to run the compendia jobs.
At the time of writing, compendia jobs require 180GB of RAM and m5.12xlarge has 192GM of RAM so it is sufficiently large to run the jobs.


### Running Tximport Early

Normally we wait until ever sample in an experiment has had Salmon run on it before we run Tximport.
However Salmon won't work on every sample, so some experiments are doomed to never make it to 100% completion.
Tximport can be run on such an experiment with:

```bash
nomad job dispatch -meta EXPERIMENT_ACCESSION=SRP009841 TXIMPORT
```

Note that if the experiment does not have at least 25 samples with at least 80% of them processed, this will do nothing.


### Checking on Local Jobs

_Note:_ The following instructions assume you have set the environment
variable NOMAD_ADDR to include the IP address of your development
machine. This can be done with:

```bash
source ./scripts/common.sh && export NOMAD_ADDR=http://$(get_ip_address):4646
```

To check on the status of a job, run:

```bash
nomad status
```

It should output something like:

```
ID                                       Type                 Priority  Status   Submit Date
DOWNLOADER                               batch/parameterized  50        running  01/31/18 18:34:05 EST
DOWNLOADER/dispatch-1517441663-4b02e7a3  batch                50        dead     01/31/18 18:34:23 EST
PROCESSOR                                batch/parameterized  50        running  01/31/18 18:34:05 EST
```

The rows whose `ID`s are `DOWNLOADER` or `PROCESSOR` are the parameterized
jobs which are waiting to dispatch Refine.bio jobs. If you don't understand
what that means, don't worry about it. All you really need to do is select
one of the jobs whose ID contains `dispatch` and whose `Submit Date`
matches the time when the job you want to check on was run, copy that full ID
(in this case `DOWNLOADER/dispatch-1517437920-ae8b77a4`), and paste it
after the previous command, like so:

```bash
nomad status DOWNLOADER/dispatch-1517441663-4b02e7a3
```

This will output a lot of information about that `Nomad Dispatch Job`,
of which we're mostly interested in the section titled **Allocations**.
Here is an example:

```
Allocations
ID        Node ID   Task Group  Version  Desired  Status    Created At
b30e4edd  fda75a5a  jobs        0        run      complete  01/31/18 18:34:23 EST
```

If you paste that after the original `nomad status` command, like so:

```bash
nomad status b30e4edd
```

you'll see a lot of information about allocation, which probably isn't
what you're interested in. Instead, you should run:

```bash
nomad logs -verbose b30e4edd
```

This command will output both the stderr and stdout logs from the container
which ran that allocation. The allocation is really a Refine.bio job.


### Development Helpers

It can be useful to have an interactive Python interpreter running within the
context of the Docker container. The `scripts/run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However, it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

## Cloud Deployment

Refine.bio requires an active, credentialed AWS account with appropriate permissions to create network infrastructure, users, compute instances and databases.

Deploys are automated to run via CirlceCI whenever a signed tag starting with a `v` is pushed to either the `dev` or `master` branches (v as in version, i.e. v1.0.0).
Tags intended to trigger a staging deploy MUST end with `-dev`, i.e. `v1.0.0-dev`.
CircleCI runs a deploy on a dedicated AWS instance so that the Docker cache can be preserved between runs.

Instructions for setting up that instance can be found in the infrastructure/deploy_box_instance_data.sh script.

To trigger a new deploy, first see what tags already exist with `git tag --list | sort --version-sort`
We have two different version counters, one for `dev` and one for `master` so a list including things like:
* v1.1.2
* v1.1.2-dev
* v1.1.3
* v1.1.3-dev

However you may see that the `dev` counter is way ahead, because we often need more than one staging deploy to be ready for a production deploy.
This is okay, just find the latest version of the type you want to deploy and increment that to get your version.
For example, if you wanted to deploy to staging and the above versions were the largest that `git tag --list` output, you would increment `v1.1.3-dev` to get `v1.1.4-dev`.

Once you know which version you want to deploy, say `v1.1.4-dev`, you can trigger the deploy with these commands:
```bash
git checkout dev
git pull origin dev
git tag -s v1.1.4-dev
git push origin v1.1.4-dev
```

`git tag -s v1.1.4-dev` will prompt you to write a tag message; please try to make it descriptive.

We use semantic versioning for this project so the last number should correspond to bug fixes and patches, the second middle number should correspond to minor changes that don't break backwards compatibility, and the first number should correspond to major changes that break backwards compatibility.
Please try to keep the `dev` and `master` versions in sync for major and minor versions so only the patch version gets out of sync between the two.

### Docker Images

Refine.bio uses a number of different Docker images to run different pieces of the system.
By default, refine.bio will pull images from the Dockerhub repo `ccdlstaging`.
If you would like to use images you have built and pushed to Dockerhub yourself you can pass the `-d` option to the `deploy.sh` script.

To make building and pushing your own images easier, the `scripts/update_my_docker_images.sh` has been provided.
The `-d` option will allow you to specify which repo you'd like to push to.
If the Dockerhub repo requires you to be logged in, you should do so before running the script using `docker login`.
The -v option allows you to specify the version, which will both end up on the Docker images you're building as the SYSTEM_VERSION environment variable and also will be the docker tag for the image.

`scripts/update_my_docker_images.sh` will not build the dr_affymetrix image, because this image requires a lot of resources and time to build.
It can instead be built with `./scripts/prepare_image.sh -i affymetrix -d <YOUR_DOCKERHUB_REPO>`.
WARNING: The affymetrix image installs a lot of data-as-R-packages and needs a lot of disk space to build the image.
It's not recommended to build the image with less than 60GB of free space on the disk that Docker runs on.

### Git Crypt

Secrets are stored using [git-crypt](https://www.agwa.name/projects/git-crypt/).
Team members can access secret files in the repo by running `git-crypt unlock`.

An existing team member can add a new team member who provides a GPG key `user.armor` with the following:
    1. `gpg --import user.armor`
    2. `git-crypt add-gpg-user --trusted KEYID`

Note that `git-crypt lock && git-crypt unlock` will reset permission secret files.
For the ssh key, this will require running `chmod 600 infrastructure/data-refinery-key.pem` before sshing onto AWS instances.

If you are adding a member to git-crypt, you should also consider adding thier GPG key to the `keys/` directory.
This directory is used to validate signed tags so that we know only trusted members are pushing deploys.


### Terraform

Once you have Terraform installed, `git-crypt` unlocked, and your AWS account credentials installed, you're ready to deploy. The correct way to deploy to the cloud is by running the `deploy.sh` script. This script will perform additional
configuration steps, such as setting environment variables, setting up Nomad job specifications, and performing database migrations. It can be used from the `infrastructure` directory like so:

```bash
./deploy.sh -u myusername -e dev -r us-east-1 -v v1.0.0 -d my-dockerhub-repo
```

This will spin up the whole system. It will usually take about 15 minutes, most of which is spent waiting for the Postgres instance to start.
The command above would spin up a development stack in the `us-east-1` region where all the resources' names would end with `-myusername-dev`.
All of the images used in that stack would come from `my-dockerhub-repo` and would be tagged with `v1.0.0`.

The `-e` specifies the environment you would like to spin up. You may specify, `dev`, `staging`, or `prod`. `dev` is meant for individuals to test infrastructure changes or to run large tests. `staging` is to test the overall system before re-deploying to `prod`.


To see what's been created at any time, you can:
```
terraform state list
```

If you want to change a single entity in the state, you can use

```
terraform taint <your-entity-from-state-list>
```

And then rerun `deploy.sh` with the same parameters you originally ran it with.


### Running Jobs

Jobs can be submitted via Nomad, either from a server/client or a local machine if you supply a server address and have an open network ingress.

To start a job with a file located on the foreman docker image:

```
nomad job dispatch -meta FILE=NEUROBLASTOMA.txt SURVEYOR_DISPATCHER
```

or to start a job with a file located in S3:

```
nomad job dispatch -meta FILE=s3://data-refinery-test-assets/NEUROBLASTOMA.txt SURVEYOR_DISPATCHER
```

### Log Consumption

All of the different Refine.bio subservices log to the same AWS CloudWatch Log
Group. If you want to consume these logs, you can use the `awslogs` tool, which
can be installed from `pip` like so:

```bash
pip install awslogs
```

or, for OSX El Capitan:

```bash
pip install awslogs --ignore-installed six
```

Once `awslogs` is installed, you can find your log group with:

```bash
awslogs groups
```

Then, to see all of the logs in that group for the past day, watching as they come in:

```bash
awslogs get <your-log-group> ALL --start='1 days' --watch
```

You can also apply a filter on these logs like so:

```bash
awslogs get <your-log-group> ALL --start='1 days' --watch --filter-pattern="DEBUG"
```

Or, look at a named log stream (with or without a wildcard.) For instance:
(Unfortunately this feature seems to be broken at the moment: https://github.com/jorgebastida/awslogs/issues/158)

```bash
awslogs get data-refinery-log-group-myusername-dev log-stream-api-nginx-access-* --watch
```

will show all of the API access logs made by Nginx.

### Dumping and Restoring Database Backups

Automatic snapshots are created automatically by RDS. Manual database dumps can be created by priveledged users with [these instructions](https://gist.github.com/syafiqfaiz/5273cd41df6f08fdedeb96e12af70e3b). Postgres versions on the host (I suggest the PGBouncer instance) must match the RDS instance version:

```bash
sudo add-apt-repository "deb http://apt.postgresql.org/pub/repos/apt/ $(lsb_release -sc)-pgdg main"
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
sudo apt-get update
sudo apt-get install postgresql-9.6
```

Archival dumps can also be provided upon request.

Dumps can be restored locally by copying the `backup.sql` file to the `volumes_postgres` directory, then executing:

```bash
docker exec -it drdb /bin/bash
psql --user postgres -d data_refinery -f /var/lib/postgresql/data/backup.sql
```

This can take a long time (>30 minutes)!

### Tearing Down

A stack that has been spun up via `deploy.sh -u myusername -e dev` can be taken down with `destroy_terraform.sh  -u myusername -e dev`.
The same username and environment must be passed into `destroy_terraform.sh` as were used to run `deploy.sh` either via the -e and -u options or by specifying `TF_VAR_stage` or `TF_VAR_user` so that the script knows which to take down.
Note that this will prompt you for confirmation before actually destroying all of your cloud resources.

## Support

Refine.bio is supported by
[Alex's Lemonade Stand Foundation](https://www.alexslemonade.org/),
with some initial development supported by the Gordon and Betty Moore
Foundation via GBMF 4552 to Casey Greene.

## Meta-README

The table of contents for this README is generated using `doctoc`.
`doctoc` can be installed with: `sudo npm install -g doctoc`
Once `doctoc` is installed the table of contents can be re-generated with: `doctoc README.md`

## License

BSD 3-Clause License.

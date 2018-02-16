# Data Refinery [![Build Status](https://circleci.com/gh/data-refinery/data-refinery/tree/dev.svg?&style=shield)](https://circleci.com/gh/data-refinery/data-refinery/)

<!-- This section needs to be drastically improved -->
Data Refinery harmonizes petabytes of publicly available biological data into
ready-to-use datasets for cancer researchers and AI/ML scientists.

The Data Refinery currently has four sub-projects contained within this repo:
- [common](./common) Contains code needed by both `foreman` and `workers`.
- [foreman](./foreman) Discovers data to download/process and manages jobs.
- [workers](./workers) Runs Downloader and Processor jobs.
- [terraform](./terraform) Manages infrastructure for the Data Refinery.

## 1. Development

### 1.1 Git Workflow

`data-refinery` uses a
[feature branch](http://nvie.com/posts/a-successful-git-branching-model/)
based workflow. New features should be developed on new feature branches, and
pull requests should be sent to the `dev` branch for code review. Merges into
`master` happen at the end of sprints, and tags in `master` correspond to
production releases.

### 1.2 Installation

To run the Data Refinery locally, you will need to have the
prerequisites installed onto your local machine. This will vary depending on
whether you are developing on a Mac or a Linux machine. Linux instructions
have been tested on Ubuntu 16.04 or later, but other Linux distributions
_should_ be able to run the necessary services. Microsoft Windows is currently
unsupported by this project.

#### 1.2.1 Linux

The following services will need to be installed:
- [Python3 and Pip]: `sudo apt-get -y install python3-pip`
- [Docker](https://www.docker.com/community-edition): Be sure to follow the
[post installation steps]
(https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user)
so Docker does not need sudo permissions.
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- git-crypt

Instructions for installing Docker and Nomad can be found by
following the link for each service. git-crypt can be installed via
`sudo apt-get install git-crypt`.

#### 1.2.2 Mac

The following services will need to be installed:
- [Docker for Mac](https://www.docker.com/docker-mac)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- [Homebrew](https://brew.sh/)
- git-crypt
- iproute2mac

Instructions for installing Docker, Nomad, and Homebrew can be found by
following the link for those services. The last three on that list can
be installed by running: `brew install iproute2mac git-crypt`.

#### 1.2.3 Virtual Environment

Run `./create_virtualenv.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire `data_refinery`
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

### 1.3 Common Dependecies

The [common](./common) sub-project contains common code which is
depended upon by the other sub-projects. So before anything else you
should prepare the distribution directory `common/dist` with this
command:

```bash
(cd common && python setup.py sdist)
```

### 1.4 Services

`data-refinery` also depends on Postgres and Nomad. Postgres can be
run in a local Docker container, but Nomad must be run on your
development machine.

#### 1.4.1 Nomad

Similarly, you will need to run a local
[Nomad](https://www.nomadproject.io/) service in development
mode. Assuming you have followed the [installation instructions](#installation), you
can do so with:

```bash
sudo -E ./run_nomad.sh
```

(_Note:_ This step may take some time because it downloads lots of files.)

Nomad is an orchestration tool which the Data Refinery uses to run
`Downloader` and `Processor` jobs. Jobs are queued by sending a message to
the Nomad agent, which will then launch a Docker container which runs
the job.

#### 1.4.2 Postgres

To start a local Postgres server in a Docker container, use:

```bash
./run_postgres.sh
```

Then, to initialize the database, run:

```bash
./common/install_db_docker.sh
```

Finally, to make the migrations to the database, use:

```bash
./common/make_migrations.sh
```

If you need to access a `psql` shell for inspecting the database, you can use:

```bash
./run_psql_shell.sh
```

### 1.5 Running Locally

Once you've built the `common/dist` directory and have
the Nomad and Postgres services running, you're ready to run
jobs. There are three kinds of jobs within the Data Refinery.

#### 1.5.1 Surveyor Jobs

Surveyor Jobs discover samples to download/process along with
recording metadata about the samples. A Surveyor Job should queue
`Downloader Jobs` to download the data it discovers.

The Surveyor can be run with the `./foreman/run_surveyor.sh`
script. The first argument to this script is the type of Surveyor Job
to run. The three valid options are:
- `survey_array_express`
- `survey_sra`
- `survey_transcriptome`

Each Surveyor Job type expects unique arguments. Details on these
arguments can be viewed by running:

```bash
./foreman/run_surveyor.sh <JOB_TYPE> -h
```

Templates and examples of valid commands to run the different types of
Surveyor Jobs are:

(1) The [Array Express](https://www.ebi.ac.uk/arrayexpress/) Surveyor
expects a single accession code:

```bash
./foreman/run_surveyor.sh survey_array_express <ARRAY_EXPRESS_ACCESSION_CODE>
```

Example:
```bash
./foreman/run_surveyor.sh survey_array_express E-MTAB-3050
```

(2) The [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) Surveyor expects a
range of SRA accession codes:

```bash
./foreman/run_surveyor.sh survey_sra <START> <END>
```

Example:
```bash
./foreman/run_surveyor.sh survey_sra DRR002116 DRR002116
```

(3) The Index Refinery Surveyor expects an
[Ensembl](http://ensemblgenomes.org/) divsion and a number of
organisms to survey:

```bash
./foreman/run_surveyor.sh survey_transcriptome <DIVISION> <NUMBER_OF_ORGANISMS>
```

Example:
```bash
./foreman/run_surveyor.sh survey_transcriptome Ensembl 1
```

#### 1.5.2 Downloader Jobs

Downloader Jobs will be queued automatically when `Surveyor Jobs`
discover new samples. However, if you just want to queue a `Downloader Job`
without running the `Surveyor`, the following command will queue a
`Downloader Job` which will download a sample from Array Express:

```bash
./workers/tester.sh queue_downloader
```

#### 1.5.3 Processor Jobs

Processor Jobs will be queued automatically by successful `Downloader Jobs`.
However, if you just want to run a `Processor Job` without first running
a `Downloader Job`, the following command will do so:

```bash
./workers/tester.sh queue_processor <PROCESSOR_TYPE>
```

Examples:
```bash
./workers/tester.sh queue_processor SRA
./workers/tester.sh queue_processor TRANSCRIPTOME_INDEX
```

#### 1.5.4 Checking on Local Jobs

_Note:_ The following instructions assume you have set the
environment variable $HOST_IP to the IP address of your
development machine. This can be done with:

```bash
source common.sh && export HOST_IP=$(get_ip_address)
```

To check on the status of a job, run:

```bash
nomad status -address http://$HOST_IP:4646
```

It should output something like:

```
ID                                       Type                 Priority  Status   Submit Date
DOWNLOADER                               batch/parameterized  50        running  01/31/18 18:34:05 EST
DOWNLOADER/dispatch-1517441663-4b02e7a3  batch                50        dead     01/31/18 18:34:23 EST
PROCESSOR                                batch/parameterized  50        running  01/31/18 18:34:05 EST
```

The rows whose `ID`s are `DOWNLOADER` or `PROCESSOR` are the parameterized
jobs which are waiting to dispatch Data Refinery jobs. If you don't understand
what that means, don't worry about it. All you really need to do is select
one of the jobs whose ID contains `dispatch` and whose `Submit Date`
matches the time when the job you want to check on was run, copy that full ID
(in this case `DOWNLOADER/dispatch-1517437920-ae8b77a4`), and paste it
after the previous command, like so:

```bash
nomad status -address http://$HOST_IP:4646 DOWNLOADER/dispatch-1517441663-4b02e7a3
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
nomad status -address http://$HOST_IP:4646 b30e4edd
```

you'll see a lot of information about allocation, which probably isn't
what you're interested in. Instead, you should run:

```bash
nomad logs -verbose -address http://$HOST_IP:4646 b30e4edd
```

This command will output both the stderr and stdout logs from the container
which ran that allocation. The allocation is really a Data Refinery job.

### 1.6 Testing

To run the entire test suite:

```bash
./run_all_tests.sh
```

These tests will also be run continuosly for each commit via CircleCI.

### 1.7 Production Deployment

_TODO_

### 1.8 Development Helpers

It can be useful to have an interactive Python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However, it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

### 1.9 Style

R files in this repo follow
[Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml).
Python Files in this repo follow
[PEP 8](https://www.python.org/dev/peps/pep-0008/). All files (including
Python and R) have a line length limit of 100 characters.

A `setup.cfg` file has been included in the root of this repo which specifies
the line length limit for the autopep8 and flake8 linters. If you run either
linter within the project's directory tree, it will enforce a line length limit
of 100 instead of 80. This will also be true for editors which rely on either
linter.

## 2. Support

`data-refinery` is supported by
[Alex's Lemonade Stand Foundation](https://www.alexslemonade.org/),
with some initial development supported by the Gordon and Betty Moore
Foundation via GBMF 4552 to Casey Greene.

## 3. License

BSD 3-Clause License.

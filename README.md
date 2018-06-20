# Refine.bio [![Build Status](https://circleci.com/gh/AlexsLemonade/refinebio/tree/dev.svg?&style=shield)](https://circleci.com/gh/AlexsLemonade/refinebio/)

<!-- This section needs to be drastically improved -->
Refine.bio harmonizes petabytes of publicly available biological data into
ready-to-use datasets for cancer researchers and AI/ML scientists.

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
    - [Linux](#linux)
    - [Mac](#mac)
    - [Virtual Environment](#virtual-environment)
    - [Common Dependecies](#common-dependecies)
    - [Services](#services)
      - [Nomad](#nomad)
      - [Postgres](#postgres)
  - [Testing](#testing)
  - [Development Helpers](#development-helpers)
  - [Style](#style)
  - [Gotchas](#gotchas)
- [Running Locally](#running-locally)
  - [Surveyor Jobs](#surveyor-jobs)
    - [ArrayExpress](#arrayexpress)
    - [Sequence Read Archive (SRA)](#sequence-read-archive-sra)
    - [Gene Expression Omnibus (GEO)](#gene-expression-omnibus-geo)
    - [Ensembl Indexes](#ensembl-indexes)
  - [Downloader Jobs](#downloader-jobs)
  - [Processor Jobs](#processor-jobs)
  - [Checking on Local Jobs](#checking-on-local-jobs)
  - [Testing](#testing-1)
    - [API](#api)
    - [Common](#common)
    - [Foreman](#foreman)
    - [Workers](#workers)
  - [Development Helpers](#development-helpers-1)
  - [Style](#style-1)
- [Production Deployment](#production-deployment)
  - [Terraform](#terraform)
  - [Autoscaling and Setting Spot Prices](#autoscaling-and-setting-spot-prices)
  - [Running Jobs](#running-jobs)
  - [Log Consumption](#log-consumption)
- [Support](#support)
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

#### Linux

The following services will need to be installed:
- [Python3 and Pip]: `sudo apt-get -y install python3-pip`
- [Docker](https://www.docker.com/community-edition): Be sure to follow the
[post installation steps]
(https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user)
so Docker does not need sudo permissions.
- [Terraform](https://www.terraform.io/)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries) can be installed on Linux clients with `sudo ./install_nomad.sh`.
- git-crypt
- jq
- iproute2

Instructions for installing Docker and Nomad can be found by
following the link for each service. git-crypt, jq, and iproute2 can be installed via
`sudo apt-get install git-crypt jq iproute2`.

When installing pip packages later in the install, you might get an error saying you need sudo permissions.
In order to fix this you have to edit your `~/.config/pip/pip.conf` to add this:

```
[install]
user = yes
no-binary = :all:
```

This sets pip to install all packages in your user directory so sudo is not required for pip intalls.

#### Mac

The following services will need to be installed:
- [Docker for Mac](https://www.docker.com/docker-mac)
- [Terraform](https://www.terraform.io/)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- [Homebrew](https://brew.sh/)
- git-crypt
- iproute2mac
- jq

Instructions for installing Docker, Nomad, and Homebrew can be found by
following the link for those services. The others on that list can
be installed by running: `brew install iproute2mac git-crypt terraform jq`.

Many of the computational processes running are very memory intensive. You will need
to [raise the amount of virtual memory available to
Docker](https://docs.docker.com/docker-for-mac/#advanced) from the default of
2GB to 12GB or 24GB, if possible.

#### Virtual Environment

Run `./create_virtualenv.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire `refinebio`
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

#### Common Dependecies

The [common](./common) sub-project contains common code which is
depended upon by the other sub-projects. So before anything else you
should prepare the distribution directory `common/dist` with this
command:

```bash
(cd common && python setup.py sdist)
```

#### Services

`refinebio` also depends on Postgres and Nomad. Postgres can be
run in a local Docker container, but Nomad must be run on your
development machine.

##### Nomad

Similarly, you will need to run a local
[Nomad](https://www.nomadproject.io/) service in development
mode. Assuming you have followed the [installation instructions](#installation), you
can do so with:

```bash
sudo -E ./run_nomad.sh
```

(_Note:_ This step may take some time because it downloads lots of files.)

Nomad is an orchestration tool which Refine.bio uses to run
`Downloader` and `Processor` jobs. Jobs are queued by sending a message to
the Nomad agent, which will then launch a Docker container which runs
the job. If address conflicts emerge, old Docker containers can be purged
with `docker container prune -f`.

##### Postgres

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

Note: there is a small chance this might fail with a `can't stat`, error. If this happens, you have
to manually change permissions on the volumes directory with `sudo chmod -R 740 volumes_postgres`
then re-run the migrations.

If you need to access a `psql` shell for inspecting the database, you can use:

```bash
./run_psql_shell.sh
```

### Testing

The end to end tests require a separate Nomad client to be running so
that the tests can be run without interfering with local
development. The second Nomad client can be started with:

```bash
sudo -E ./run_nomad.sh -e test
```

To run the entire test suite:

```bash
./run_all_tests.sh
```

(_Note:_ Running all the tests can take some time because it downloads a lot of files)

These tests will also be run continuosly for each commit via CircleCI.

### Development Helpers

It can be useful to have an interactive Python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However, it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

### Style

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

### Gotchas

During development, you make encounter some occasional strangeness. Here's
some things to watch out for:

  - If builds are failing, increase the size of Docker's memory allocation.
  - If Docker images are failing mysteriously during creation, it may
be the result of Docker's `Docker.qcow2` or `Docker.raw` file filling. You
can prune old images with `docker system prune -a`.
  - If it's killed abruptly, the containerized Postgres images can be
  left in an unrecoverable state. Annoying.
  - Since we use multiple Docker instances, don't forget to `./update_models`

## Running Locally

Once you've built the `common/dist` directory and have
the Nomad and Postgres services running, you're ready to run
jobs. There are three kinds of jobs within Refine.bio.

### Surveyor Jobs

Surveyor Jobs discover samples to download/process along with
recording metadata about the samples. A Surveyor Job should queue
`Downloader Jobs` to download the data it discovers.

The Surveyor can be run with the `./foreman/run_surveyor.sh`
script. The first argument to this script is the type of Surveyor Job
to run. The valid options are:
- `survey_all`
- `survey_array_express`
- `survey_sra`
- `survey_geo`
- `survey_transcriptome`

Each Surveyor Job type expects unique arguments. Details on these
arguments can be viewed by running:

```bash
./foreman/run_surveyor.sh <JOB_TYPE> -h
```

You can also supply a newline-deliminated file to `survey_all` which will
dispatch survey jobs based on acession codes like so:

Example:
```bash
./foreman/run_surveyor.sh survey_all --file MY_BIG_LIST_OF_CODES.txt
```

Templates and examples of valid commands to run the different types of
Surveyor Jobs are:

#### ArrayExpress

The [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) Surveyor
expects a single accession code:

```bash
./foreman/run_surveyor.sh survey_array_express --accession <ARRAY_EXPRESS_ACCESSION_CODE>
```

Example:
```bash
./foreman/run_surveyor.sh survey_array_express --accession E-MTAB-3050
```

You can also supply a list in a file:
```bash
./foreman/run_surveyor.sh survey_array_express --file my_neuroblastoma_list.txt
```

#### Sequence Read Archive (SRA)

The [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) Surveyor expects a
range of SRA accession codes:

```bash
./foreman/run_surveyor.sh survey_sra --accession <ACCESSION>
```

Example (single read):
```bash
./foreman/run_surveyor.sh survey_sra --accession DRR002116
```

Example (paired read):
```bash
./foreman/run_surveyor.sh survey_sra --accession SRR6718414
```

#### Gene Expression Omnibus (GEO)

The GEO surveyor expects an exession code or a file:

```bash
./foreman/run_surveyor.sh survey_geo --accession <ARRAY_EXPRESS_ACCESSION_CODE> --file <FILEPATH.txt>
```

Example:
```bash
./foreman/run_surveyor.sh survey_geo --file NEUROBLASTOMA.txt
```

#### Ensembl Indexes

The Index Refinery Surveyor expects an [Ensembl](http://ensemblgenomes.org/) divsion and a number of
organisms to survey:

```bash
./foreman/run_surveyor.sh survey_transcriptome <DIVISION> <NUMBER_OF_ORGANISMS>
```

Example:
```bash
./foreman/run_surveyor.sh survey_transcriptome Ensembl 1
```

A specific organism can also be specified:

```bash
./foreman/run_surveyor.sh survey_transcriptome Ensembl 2 "Homo Sapiens"
```

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

Or for more information run:
```bash
./workers/tester.sh -h
```

### Checking on Local Jobs

_Note:_ The following instructions assume you have set the environment
variable NOMAD_ADDR to include the IP address of your development
machine. This can be done with:

```bash
source common.sh && export NOMAD_ADDR=http://$(get_ip_address):4646
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

### Testing

To run the entire test suite:

```bash
./run_all_tests.sh
```

These tests will also be run continuosly for each commit via CircleCI.

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


### Development Helpers

It can be useful to have an interactive Python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However, it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

### Style

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

## Production Deployment

Refine.bio requires an active, credentialed AWS account with appropriate permissions to create network infrastructure, users, compute instances and databases.

### Terraform

Once you have Terraform installed and your AWS account credentials installed, you can plan your terraform deployment like so (from the `infrastructure` directory):

```bash
TF_VAR_user=myusername TF_VAR_stage=dev TF_VAR_region=us-east-1 terraform plan
```

If that worked fine, then to deploy:

```bash
TF_VAR_user=myusername TF_VAR_stage=dev TF_VAR_region=us-east-1 ./deploy.sh
```

This will spin up the whole system. It will usually take about 15 minutes, most of which is spent waiting for the Postgres instance to start.

To see what's been created at any time, you can:
```
terraform state list
```

If you want to change a single entity in the state, you can use

```
terraform taint <your-entity-from-state-list>; tf plan; tf apply;
```

To tear down the entire system:

```
terraform destroy
```

For convenience, a `deploy.sh` script is also provided, which will perform additional
steps to configure (such as setting up Nomad job specifications and performing database migrations) and prepare the entire system. It can be used simply (from the `infrastructure` directory), like so:

```
./deploy.sh
```

### Autoscaling and Setting Spot Prices

`refinebio` uses AWS Auto Scaling Groups to provide elastic capacity for large
work loads. To do this, we use "Spot Requests". To find a good bid price for
your instance type, use the [spot request
page](https://console.aws.amazon.com/ec2sp/v1/spot/home?region=us-east-1) and
then click on the **view pricing history** chart. Choose your instance type and
then choose a bid price that is slightly higher than the current price for your
availability zone (AZ). [This graph](https://cdn-images-1.medium.com/max/1600/0*gk64fOrhSFBoFGXK.) is useful for understanding instance types.

Then set your `TF_VAR_client_instance_type`, `TF_VAR_spot_price` and
`TF_VAR_max_clients` to configure your scaling instance types, cost and size.
`TF_VAR_scale_up_threshold` and `TF_VAR_scale_down_threshold` define the queue
lengths which trigger the scaling alarms, though you probably won't need to
tweak these as much.

### Running Jobs

Jobs can be submitted via Nomad, either from a server/client or a local machine if you supply a server address and have an open network ingress.

To start the Neuroblastoma job:

```
nomad job dispatch -meta COMMAND=survey_array_express -meta FILE=NEUROBLASTOMA.txt SURVEYOR
```

or the Zebrafish job:

```
nomad job dispatch -meta COMMAND=survey_all -meta FILE=s3://data-refinery-test-assets/ZEBRAFISH.txt SURVEYOR
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

```bash
awslogs get data-refinery-log-group-myusername-dev log-stream-api-nginx-access-* --watch
```

will show all of the API access logs made by Nginx.

## Support

Refine.bio is supported by
[Alex's Lemonade Stand Foundation](https://www.alexslemonade.org/),
with some initial development supported by the Gordon and Betty Moore
Foundation via GBMF 4552 to Casey Greene.

## License

BSD 3-Clause License.

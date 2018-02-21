# Data Refinery [![Build Status](https://circleci.com/gh/data-refinery/data-refinery/tree/dev.svg?&style=shield)](https://circleci.com/gh/data-refinery/data-refinery/)

<!-- This section needs to be drastically improved -->
Data Refinery harmonizes petabytes of publicly available biological data into ready-to-use datasets for cancer researchers and AI/ML scientists.

The Data Refinery currently has four sub-projects contained within this repo:
- [common](./common) Contains code needed by both `foreman` and `workers`.
- [foreman](./foreman) Discovers data to download/process and manages jobs.
- [workers](./workers) Runs Downloader and Processor jobs.
- [terraform](./terraform) Manages infrastructure for the Data Refinery.

## Development

### Git Workflow

`data-refinery` uses a [feature branch](http://nvie.com/posts/a-successful-git-branching-model/) based workflow. New features should be
developed on new feature branches, and pull requests should be sent to
the `dev` branch for code review. Merges into `master` happen at the end
of sprints, and tags in `master` correspond to production releases.


### Installation

To run the Data Refinery locally you will need to have the
prerequisites installed onto your local machine. This will vary
whether you are developing on a Mac or a Linux machine. Linux
instructions have been tested on Ubuntu, but other Linux distributions
_should_ be able to run the necessary services. Windows is currently
unsupported by this project.

#### Linux

The following services will need to be installed:
- [Python3 and Pip]: `sudo apt-get -y install python3-pip`
- [Docker](https://www.docker.com/community-edition): Be sure to follow the [post installation steps](https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user) so Docker does not need sudo permissions.
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- git-crypt

Instructions for installing Docker and Nomad can be found by
following the link for each service. git-crypt can be installed via
`sudo apt-get install git-crypt`.

#### Mac
The following services will need to be installed:
- [Docker for Mac](https://www.docker.com/docker-mac)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- [Homebrew](https://brew.sh/)
- git-crypt
- iproute2mac

Instructions for installing Docker, Nomad, and Homebrew can be found by
following the link for those services. The last three on that list can
be installed by running: `brew install iproute2mac git-crypt`.

Many of the computational processes running are very memory intensive. You will need
to [raise the amount of virtual memory available to Docker](https://docs.docker.com/docker-for-mac/#advanced) from the default of 2GB to 12GB or 24GB, if possible. 

#### Virtual Environment

Run `./create_virtualenv.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire `data_refinery`
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

### Common Dependecies

The [common](./common) sub-project contains common code which is
depended upon by the other sub-projects. So before anything else you
should prepare the distribution directory for `common` with this
command:

```bash
(cd common && python setup.py sdist)
```

### Services

`data-refinery` also depends on Postgres and Nomad. Postgres can be
run in a local Docker container, but Nomad must be run on your
development machine.


#### Nomad

Similarly, you will need to run a local
[Nomad](https://www.nomadproject.io/) service in development
mode. Assuming you have followed the [installation instructions](#installation), you
can do so with:

```bash
sudo -E ./run_nomad.sh
```

Nomad is an orchestration tool which the Data Refinery uses to run
Downloader and Processor jobs. Jobs are queued by sending a message to
the Nomad agent, which will then launch a Docker container which runs
the job. If address conflicts emerge, old Docker containers can be purged
with `docker container prune -f`.


#### Postgres

To start a local Postgres server in a Docker container, use:

```bash
./run_postgres.sh
```

Then, to initalize the database, run:

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

### Running Locally

Once you've run built the `common` distribution directory and you have
the Nomad and Postgres services running, you're ready to run
jobs. There are three kinds of jobs within the Data Refinery.

#### Surveyor Jobs

Surveyor Jobs discover samples to download/process along with
recording metadata about the samples. A Surveyor Job should queue
Downloader Jobs to download the data it discovers.

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

The [Array Express](https://www.ebi.ac.uk/arrayexpress/) Surveyor
expects a single accession code:

```bash
./foreman/run_surveyor.sh survey_array_express <ARRAY_EXPRESS_ACCESSION_CODE>
```

Example:
```bash
./foreman/run_surveyor.sh survey_array_express E-MTAB-3050
```

The [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) Surveyor expects a
range of SRA accession codes:

```bash
./foreman/run_surveyor.sh survey_sra <START> <END>
```

Example (single read):
```bash
./foreman/run_surveyor.sh survey_sra DRR002116
```

Example (paired read):
```bash
./foreman/run_surveyor.sh survey_sra SRR6718414
```

The Index Refinery Surveyor expects an
[Ensembl](http://ensemblgenomes.org/) divsion and a number of
organisms to survey:

```bash
./foreman/run_surveyor.sh survey_transcriptome <DIVISION> <NUMBER_OF_ORGANISMS>
```

Example:
```bash
./foreman/run_surveyor.sh survey_transcriptome Ensembl 1
```

A specific organism can also be specified:

Example:
```bash
./foreman/run_surveyor.sh survey_transcriptome Ensembl 2 "Homo Sapiens"
```

#### Downloaders

Downloader Jobs will be queued automatically when Surveyor Jobs
discover new samples. However if you just want to queue a Downloader
Job without running the Surveyor, the following command will queue a
Downloader Job which will download a sample from Array Express:

```bash
./workers/tester.sh queue_downloader
```

For example:
```bash
./common/make_migrations.sh; cd common; python setup.py sdist; cd ..; ./workers/tester.sh run_downloader_job --job-name=ARRAY_EXPRESS --job-id=1
```


#### Processors

Processor Jobs will be queued automatically by successful Downloader
Jobs. However, if you just want to run a Processor Job without first
needing to run a Downloader Job, the following command will do so:

```bash
./workers/tester.sh queue_processor <PROCESSOR_TYPE>
```

An individual processor job can be run with:
```bash
./workers/tester.sh queue_processor SRA
./workers/tester.sh queue_processor TRANSCRIPTOME_INDEX
```


#### Checking on Local Jobs

_Note:_ The following instructions assume you have set the
environment variable $HOST_IP to the the IP address of your
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

The DOWNLOADER and PROCESSOR IDs are the parameterized jobs which are
waiting to dispatch Data Refinery jobs. If you don't understand what
that means, don't worry about it. All you really need to do is select
one of the IDs which contains `dispatch` and has a `Submit Date`
matching when the job you want to check on was ran. Copy that full ID,
in this case `DOWNLOADER/dispatch-1517437920-ae8b77a4`, and paste it
after the previous command like so:

```bash
nomad status -address http://$HOST_IP:4646 DOWNLOADER/dispatch-1517441663-4b02e7a3
```

This will output a lot of information about that Nomad Dispatch Job,
of which we're mostly interested in the section titled
**Allocations**, which looks like this:

```
Allocations
ID        Node ID   Task Group  Version  Desired  Status    Created At
b30e4edd  fda75a5a  jobs        0        run      complete  01/31/18 18:34:23 EST
```

If you paste that after the original `nomad status` command like so:

```bash
nomad status -address http://$HOST_IP:4646 b30e4edd
```

You'll see a lot of information about allocation, which probably isn't
what you're interested in. Instead you should run:

```bash
nomad logs -verbose -address http://$HOST_IP:4646 b30e4edd
```

Which will output both the stderr and stdout logs from the container
which ran that allocation, which is really a Data Refinery job.


### Testing

To run the test entire suite:

```bash
./run_all_tests.sh
```

These tests will also be run continuosly for each commit via CircleCI.

### Production Deployment

_TODO_


### Development Helpers

It can be useful to have an interactive python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

### Style

R files in this repo follow
[Google's R Style Guide](https://google.github.io/styleguide/Rguide.xml).
Python Files in this repo follow
[PEP 8](https://www.python.org/dev/peps/pep-0008/). All files (including
python and R) have a line limit of 100 characters.

A `setup.cfg` file has been included in the root of this repo which specifies
the line length limit for the autopep8 and flake8 linters. If you run either
of those programs from anywhere within the project's directory tree they will
enforce a limit of 100 instead of 80. This will also be true for editors which
rely on them.


## Support

`data-refinery` is supported by [Alex's Lemonade Stand
Foundation](https://www.alexslemonade.org/), with some initial
development supported by the Gordon and Betty Moore Foundation via
GBMF 4552 to Casey Greene.

## License

BSD 3-Clause License.

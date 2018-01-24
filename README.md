# Data Refinery [![Build Status](ttps://circleci.com/gh/data-refinery/data-refinery/tree/dev.svg?&style=shield)](https://circleci.com/gh/data-refinery/data-refinery/)

A project to download, process, aggregate, and serve bioinformatic data now
supported by Alex's Lemonade Stand Foundation, with some initial development
supported by the Gordon and Betty Moore Foundation via GBMF 4552 to Casey
Greene.

## Development

### Git Workflow

`data-refinery` uses a [feature branch](http://nvie.com/posts/a-successful-git-branching-model/ based workflow. New features should be
developed on new feature branches, and pull requests should be sent to
the `-dev` branch for code review. Merges into `-master` happen at the end
of sprints, and tags in `master` correspond to production releases.

### Running Locally

_Note: The following steps assume you have already installed PostgreSQL (>=9.4)
and Python (>=3.5) on Ubuntu 16.04 and higher, as well as OSX. OSX users will
need [Docker for Mac](https://www.docker.com/docker-mac) installed, as well
as Homebrew, for to `brew install iproute2mac git-crypt git-lfs`._

Run `./install.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire data_refinery
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

### Services

`data-refinery` also depends on Postgres and Nomad, both of which can be run
in local Docker containers.

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

#### Nomad

Similarly, you will need to run a local Nomad service instance in a Docker
container. You can do so with:

```bash
./run_nomad.sh
```


### Testing

To run the test entire suite:

```bash
./run_all_tests.sh
```

These tests will also be run continuosly for each commit via CircleCI.

### Production Deployment

_TODO_

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

It can be useful to have an interactive python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

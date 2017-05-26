# Data Refinery

A project to download, process, aggregate, and serve bioinformatic data
supported by Greene Lab.

## Getting Started

Note: The following steps assume you have already installed PostgreSQL (>=9.4)
and Python (Most versions should work, but this has been tested with Python 3.5)
on Ubuntu (Tested with 16.04. It should be possible to use other versions or
even a Mac though).

Run `./install.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire data_refinery
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

You'll also need to set up the database. See `data_models/README.md` for
instructions on doing so.

## Development

For development it will be necessary to run a rabbitmq docker image:
```
docker run -d --hostname my-rabbit --name some-rabbit rabbitmq:3
```

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

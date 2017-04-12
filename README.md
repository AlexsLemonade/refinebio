# Data Refinery

A project to download, process, aggregate, and serve bioinformatic data
supported by Greene Lab.

## Getting Started

Note: The following steps assume you have already installed PostgreSQL (>=9.4) and Python (Most versions should work, but this has been tested with Python 3.5) on Ubuntu (Tested with 16.04. It should be possible to use other versions or even a Mac though).

Run `./install.sh` to set up the virtualenv. It will activate the `dr_env`
for you the first time. This virtualenv is valid for the entire data_refinery
repo. Sub-projects each have their own virtualenvs which are managed by their
containers. When returning to this project you should run
`source dr_env/bin/activate` to reactivate the virtualenv.

You'll also need to set up the database. See `data_models/README.md` for
instructions on doing so.

## Development

It can be useful to have an interactive python interpreter running within the
context of the Docker container. The `run_shell.sh` script has been provided
for this purpose. It is in the top level directory so that if you wish to
reference it in any integrations its location will be constant. However it
is configured by default for the Foreman project. The interpreter will
have all the environment variables, dependencies, and Django configurations
for the Foreman project. There are instructions within the script describing
how to change this to another project.

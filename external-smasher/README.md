# External Smasher
This project provides a script to allow for the processing of external data along with refine.bio
data without uploading the data to our refine.bio servers. Our servers automatically import any data
you provide us and make it available to anyone that wants it. If you want to process your data
alongside some from refine.bio but would prefer to keep your data private, you can use this
tool to process and combine it instead of going through our servers.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Installation](#installation)
  - [Linux](#linux)
  - [Mac](#mac)
  - [Common Dependecies](#common-dependecies)
  - [Services](#services)
    - [Nomad](#nomad)
    - [Postgres](#postgres)
- [Processing](#processing)
  - [Setup](#setup)
  - [Running](#running)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Installation
Our workflow relies pretty heavily on bash and docker, both things that do not play very well with
Windows. Because of that, at the moment we only support macOS and Linux.

(_Note:_ All of the following instructions expect your working directory to be the root `refinebio`
directory)

### Linux

Ubuntu 16.04 is the distro that is used most heavily by our team, and as such is the distribution we
officially support. Other distributions _should_ work, but may require more manual setup.

The following services will need to be installed:
- [Python3 and Pip]: `sudo apt-get -y install python3-pip`
- [Docker](https://www.docker.com/community-edition): Be sure to follow the
[post installation steps]
(https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user)
so Docker does not need sudo permissions.
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries) can be installed
on Linux clients with `sudo ./install_nomad.sh`.
- [jq](https://stedolan.github.io/jq/)
- [iproute2](https://wiki.linuxfoundation.org/networking/iproute2)

Instructions for installing Docker and Nomad can be found by following the link for each service.
jq, and iproute2 can be installed via `sudo apt-get install jq iproute2`.

### Mac

The following services will need to be installed:
- [Docker for Mac](https://www.docker.com/docker-mac)
- [Terraform](https://www.terraform.io/)
- [Nomad](https://www.nomadproject.io/docs/install/index.html#precompiled-binaries)
- [Homebrew](https://brew.sh/)
- [git-crypt](https://www.agwa.name/projects/git-crypt/)
- [iproute2mac](https://github.com/brona/iproute2mac)
- [jq](https://stedolan.github.io/jq/)

Instructions for installing Docker, Nomad, and Homebrew can be found by
following the link for those services. The others on that list can
be installed by running: `brew install iproute2mac git-crypt terraform jq`.

Many of the computational processes running are very memory intensive. You will need
to [raise the amount of virtual memory available to
Docker](https://docs.docker.com/docker-for-mac/#advanced) from the default of
2GB to 12GB or 24GB, if possible.

### Common Dependecies

The [common](./common) sub-project contains common code which is
depended upon by the other sub-projects. So before anything else you
should prepare the distribution directory `common/dist` with this
command:

```bash
(cd common && python setup.py sdist)
```

### Services

`refinebio` also depends on Postgres and Nomad. Postgres can be
run in a local Docker container, but Nomad must be run on your
development machine.

#### Nomad

Similarly, you will need to run a local [Nomad](https://www.nomadproject.io/) service in development
mode. Assuming you have followed the [installation instructions](#installation), you can do so with:

```bash
sudo -E ./run_nomad.sh
```

(_Note:_ This step may take some time because it downloads lots of files.)

#### Postgres

To setup the postgres database for the first time, use:
```bash
./external-smasher/setup_postgres.sh
```

After you set it up the first time, you can start it again after startup with:
```bash
./run_postgres.sh
```

Note: there is a small chance this might fail with a `can't stat`, error. If this happens, you have
to manually change permissions on the volumes directory with `sudo chmod -R 740 volumes_postgres`
then re-run the migrations.

## Processing
All of the functionality of our external processor is exposed through the shell script
`./external-smasher/combine.sh`. This script will combine data downloaded by refine.bio with a
local folder of raw data. Currently, the only raw data formats we support are Affymetrix and
RNA-SEQ, and all the samples must be from the same organism.

### Setup
Ahead of processing your local data, you must download the samples you want from refine.bio.
_Please make sure to select 'NONE' for the transformation on the data processed by refine.bio._ If
you select a different transformation, then the refine.bio data will become double-transformed
during the combination with your private data and the results __will__ be incorrect.

### Running
The script must be ran with these options:
```bash
./external-smasher/combine.sh -I REFINEBIO_DATA -D YOUR_DATA -O ORGANISM
```
`-I` specifies the path to the data downloaded from refine.bio, which must be in the `zip` format.

`-D` specifies the path to a directory holding your data. Your data must be in the Affymetrix or
RNA-SEQ formats.

`-O` specifies the name of the organism, and must be in underscore-delimited uppercase (i.e.
GALLUS_GALLUS)

There is also one optional argument:

`-T` specifies a transformation to be applied to the data after conversion. The default option is
'MINMAX', but other valid options are 'NONE', 'STANDARD', and 'ROBUST'.
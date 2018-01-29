#! /bin/bash

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# CircleCI Docker won't make this by default for some reason
if [ ! -d /tmp/volumes_postgres ]; then
  mkdir /tmp/volumes_postgres
fi
VOLUMES=/tmp/volumes_postgres

# via https://hub.docker.com/_/postgres/
# 9.6.6 is the current (as of Jan 23 2018) RDS most recent version.
# Password can be exposed to git/CI because this is only for dev/testing purposes, not real data.
docker run --name drdb -v "$VOLUMES":/var/lib/postgresql/data -e POSTGRES_PASSWORD=mysecretpassword -d postgres:9.6.6
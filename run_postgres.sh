#! /bin/bash

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd "$script_directory"

# CircleCI Docker won't make this by default for some reason
# This doubly nested directory is a hacky workaround to prevent permissions issues.
# Suggested here:
# https://github.com/docker/for-linux/issues/380#issuecomment-436419102
VOLUMES="$script_directory/volumes_postgres/volumes_postgres"
if [ ! -d "$VOLUMES" ]; then
  mkdir -p "$VOLUMES"
fi

# Check if a docker database named "drdb" exists, and if so just run it
if [[ $(docker ps -a --filter name=drdb -q) ]]; then
  docker start drdb > /dev/null
  echo "Started database."
# Otherwise, install it from docker hub
else
  # via https://hub.docker.com/_/postgres/
  # 9.6.6 is the current (as of Jan 23 2018) RDS most recent version.
  # Password can be exposed to git/CI because this is only for dev/testing purposes, not real data.
  echo "Installing database..."
  docker run -p 5432:5432 --name drdb -v "$VOLUMES":/var/lib/postgresql/data -e POSTGRES_PASSWORD=mysecretpassword -d postgres:9.6.6
fi

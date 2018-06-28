#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Ensure that Nomad is running first
# Double w in `ps` will cause the columns to never be truncated regardless of environment.
if ! ps auxww | grep test_nomad | grep -v grep > /dev/null; then
    echo "You must start the nomad test environment first with" >&2
    echo "'sudo -E ./run_nomad.sh -e test'" >&2
    exit 1
# Then ensure postgres is running
elif ! [[ $(docker ps --filter name=drdb -q) ]]; then
    echo "You must start Postgres first with './run_postgres.sh'" >&2
    exit 1
fi

./prepare_image.sh -i foreman -s foreman

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file foreman/environments/test \
       --link drdb:postgres \
       -it ccdlstaging/dr_foreman bash -c "$(run_tests_with_coverage $@)"

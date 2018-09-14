#!/bin/bash

# script for executing Django PyUnit Tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Ensure that postgres is running
if ! [[ $(docker ps --filter name=drdb -q) ]]; then
    echo "You must start Postgres first with:" >&2
    echo "./run_postgres.sh" >&2
    exit 1
fi

./prepare_image.sh -i common_tests -s common

source common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
HOST_IP=$(get_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file common/environments/test \
       --link drdb:postgres \
       -it ccdlstaging/dr_common_tests bash -c "$(run_tests_with_coverage $@)" --parallel

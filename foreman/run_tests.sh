#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# Set up the data volume directory if it does not already exist.
# Since the end-to-end tests are run from the Foreman image, use the
# top level test_volume rather than one nested within the foreman
# directory.
project_root=$(cd .. && pwd)
volume_directory="$project_root/test_volume"
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
    echo "You must start the nomad test environment first with:" >&2
    echo "sudo -E ./run_nomad.sh -e test" >&2
    exit 1
# Then ensure postgres is running
elif ! [[ $(docker ps --filter name=drdb -q) ]]; then
    echo "You must start Postgres first with:" >&2
    echo "./run_postgres.sh" >&2
    exit 1
# Then ensure elasticsearch is running
elif ! [[ $(docker ps --filter name=dres -q) ]]; then
    echo "You must start Elasticsearch first with:" >&2
    echo "./run_es.sh" >&2
    exit 1
fi

./prepare_image.sh -i foreman -s foreman

source common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
ES_HOST_IP=$(get_docker_es_ip_address)
HOST_IP=$(get_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --add-host=elasticsearch:$ES_HOST_IP \
       --env-file foreman/environments/test \
       --volume $volume_directory:/home/user/data_store \
       -it ccdlstaging/dr_foreman bash -c "$(run_tests_with_coverage $@)"

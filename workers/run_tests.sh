#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"  | xargs readlink -f`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the test data volume directory if it does not already exist
volume_directory="$script_directory/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod 775 $volume_directory
fi

docker build -t dr_worker_tests -f workers/Dockerfile.salmon_runner .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/test \
       --volume $volume_directory:/home/user/data_store \
       -i dr_worker_tests test --no-input "$@"

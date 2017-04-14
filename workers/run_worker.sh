#!/bin/bash

# Script for running a Data Refinery Worker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod 775 $volume_directory
fi

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')

# docker build -t dr_worker -f workers/Dockerfile.worker .
docker build -t dr_conda -f workers/Dockerfile.conda .
docker run \
       --link some-rabbit:rabbit \
       --name worker1 \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/dev \
       --volume $volume_directory:/home/user/data_store \
       --detach \
       dr_conda

#!/bin/bash

# Script for running the Data Refinery Surveyor container

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"  | xargs realpath`
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

docker build -t dr_foreman -f foreman/Dockerfile .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --add-host=database:$HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file foreman/environments/dev \
       --volume $volume_directory:/home/user/data_store \
       dr_foreman survey_array_express "$@"

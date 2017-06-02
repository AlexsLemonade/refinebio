#!/bin/bash

# Running this script will start an interactive python shell running
# within the context of a Docker container.
# By default the Docker container will be for the foreman project.
# This can be changed by modifying the --env-file command line arg,
# changing foreman/Dockerfile to the appropriate Dockerfile,
# changing the volume_directory path,
# and by modifying the Dockerfile.shell file appropriately.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`
cd $script_directory

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/foreman/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod 775 $volume_directory
fi

docker build -t dr_shell -f foreman/Dockerfile .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --link message-queue:rabbit \
       --add-host=database:$HOST_IP \
       --env-file foreman/environments/dev \
       --volume /tmp:/tmp \
       --volume $volume_directory:/home/user/data_store \
       --entrypoint ./manage.py \
       --interactive dr_shell shell

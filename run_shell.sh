#!/bin/bash

# Running this script will start an interactive python shell running
# within the context of a Docker container.
# By default the Docker container will be for the foreman project.
# This can be changed by modifying the --env-file command line arg
# and by modifying the Dockerfile.shell file appropriately.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"`
cd $script_directory

docker build -t dr_shell -f Dockerfile.shell .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --link some-rabbit:rabbit \
       --add-host=database:$HOST_IP \
       --env-file foreman/environments/dev \
       --volume /tmp:/tmp \
       --interactive dr_shell

#!/bin/bash

# Running this script will start an interactive python shell running
# within the context of a Docker container.
# By default the Docker container will be for the foreman project.
# This can be changed by modifying the --env-file command line arg
# and by modifying the Dockerfile.shell file appropriately.

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --link some-rabbit:rabbit \
       --add-host=database:$HOST_IP \
       --env-file /home/kurt/Development/data_refinery/foreman/environments/dev \
       -v /tmp:/tmp \
       -i dr_shell

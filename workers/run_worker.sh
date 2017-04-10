#!/bin/bash

# Script for running a Data Refinery Worker container.

# This script should always run as if it were being called from
# the directory it lives in.
cd "$( dirname "${BASH_SOURCE[0]}" )"

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t dr_worker -f workers/Dockerfile.worker .

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --link some-rabbit:rabbit \
       --name worker1 \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/dev \
       -d dr_worker

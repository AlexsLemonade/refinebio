#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t foreman_tests -f foreman/Dockerfile.tests .

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --add-host=database:$HOST_IP \
       --env-file foreman/environments/test \
       -i foreman_tests "$@"

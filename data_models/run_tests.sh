#!/bin/bash

# script for executing Django PyUnit Tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $script_directory

docker build -t models_tests -f Dockerfile.tests .

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --add-host=database:$HOST_IP \
       --env-file environments/test \
       -i models_tests "$@"

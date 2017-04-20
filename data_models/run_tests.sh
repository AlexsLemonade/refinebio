#!/bin/bash

# script for executing Django PyUnit Tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
cd `dirname "${BASH_SOURCE[0]}"`

docker build -t models_tests -f Dockerfile.tests .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')
docker run \
       --add-host=database:$HOST_IP \
       --env-file environments/test \
       -i models_tests "$@"

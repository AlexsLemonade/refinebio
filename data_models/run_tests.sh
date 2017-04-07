#!/bin/bash

# script for executing Django Unittests within Docker container.

docker build -t models_tests -f Dockerfile.tests .

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --add-host=database:$HOST_IP \
       --env-file /home/kurt/Development/data_refinery/data_models/environments/test \
       -i models_tests "$@"

#!/bin/bash

# script for executing Django Unittests within Docker container.

cd ..

docker build -t foreman_tests -f foreman/Dockerfile.tests .

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --add-host=database:$HOST_IP \
       --env-file /home/kurt/Development/data_refinery/foreman/environments/test \
       -i foreman_tests "$@"

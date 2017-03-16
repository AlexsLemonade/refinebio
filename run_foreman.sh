#!/bin/bash

# convenience script for running the foreman docker container

HOST_IP=$(ifconfig eth0 | grep "inet " | awk -F'[: ]+' '{ print $4 }')
docker run \
       --link some-rabbit:rabbit \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/dev \
       test_master

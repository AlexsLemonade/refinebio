#!/bin/bash

# Script for running a simple tester container to test the worker.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`realpath $( dirname "${BASH_SOURCE[0]}" )`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t test_master -f workers/Dockerfile .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --link some-rabbit:rabbit \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/dev \
       --entrypoint ./manage.py \
       test_master queue_test_task

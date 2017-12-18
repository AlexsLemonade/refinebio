#!/bin/bash

# Script for running a simple tester container to test the worker.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"  | xargs realpath`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t dr_tester -f workers/Dockerfile.tests .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --add-host=database:$HOST_IP \
       --env-file workers/environments/dev \
       --entrypoint ./manage.py \
       dr_worker queue_processor "TRANSCRIPTOME_INDEX"

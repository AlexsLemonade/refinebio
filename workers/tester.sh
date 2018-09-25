#!/bin/bash

# Script for running a django management command to test the worker.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..


docker build -t dr_worker -f workers/Dockerfile .

volume_directory="$script_directory/volume"

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file workers/environments/dev \
       --entrypoint ./manage.py \
       --volume $volume_directory:/home/user/data_store \
       --link drdb:postgres \
       dr_worker "$@"

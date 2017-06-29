#!/bin/bash

# Script for migrating the database using a Docker container so no
# virtual environment is needed on the host machine.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"  | xargs realpath`
echo "SCRIPT DIRECTORY!!!!!!!!!!!!!!!!!!"
echo $script_directory
cd $script_directory

docker build -t dr_models -f Dockerfile .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --volume $script_directory/data_refinery_models:/home/user/data_refinery_models \
       --add-host=database:$HOST_IP \
       --env-file environments/dev \
       --interactive \
       dr_models makemigrations data_refinery_models

docker run \
       --volume $script_directory/data_refinery_models:/home/user/data_refinery_models \
       --add-host=database:$HOST_IP \
       --env-file environments/dev \
       dr_models

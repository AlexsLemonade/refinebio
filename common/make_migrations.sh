#!/bin/bash

# Script for migrating the database using a Docker container so no
# virtual environment is needed on the host machine.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t dr_models -f common/Dockerfile .

source common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
HOST_IP=$(get_ip_address)

docker run \
       --volume $script_directory/data_refinery_common:/home/user/data_refinery_common \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file common/environments/dev \
       --interactive \
       --link drdb:postgres \
       dr_models makemigrations data_refinery_common

docker run \
       --volume $script_directory/data_refinery_common:/home/user/data_refinery_common \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file common/environments/dev \
       --link drdb:postgres \
       dr_models migrate

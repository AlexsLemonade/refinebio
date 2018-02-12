#!/bin/bash

# Script for running a Data Refinery Worker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

docker build -t dr_worker -f workers/Dockerfile .

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
       --name worker1 \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
       --env AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
       --env-file workers/environments/dev \
       --volume $volume_directory:/home/user/data_store \
       --link drdb:postgres \
       --detach \
       dr_worker

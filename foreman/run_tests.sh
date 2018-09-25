#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t dr_foreman -f foreman/Dockerfile .

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)
NOMAD_HOST_IP=$(get_docker_nomad_ip_address)
NOMAD_LINK=$(get_nomad_link_option)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$NOMAD_HOST_IP \
       --env-file foreman/environments/test \
       --link drdb:postgres $NOMAD_LINK \
       -i dr_foreman test --no-input "$@"

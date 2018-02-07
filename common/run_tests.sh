#!/bin/bash

# script for executing Django PyUnit Tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t common_tests -f common/Dockerfile .

source common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
NOMAD_HOST_IP=$(get_docker_nomad_ip_address)
HOST_IP=$(get_ip_address)
NOMAD_LINK=$(get_nomad_link_option)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$NOMAD_HOST_IP \
       --env-file common/environments/test \
       --link drdb:postgres $NOMAD_LINK \
       -i common_tests bash -c 'coverage run --source="." manage.py test --no-input "$@"; coverage report -m'

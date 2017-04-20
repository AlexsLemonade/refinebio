#!/bin/bash

# Script for running the Data Refinery Surveyor container

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

docker build -t dr_surveyor -f foreman/Dockerfile.surveyor .

HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

docker run \
       --link some-rabbit:rabbit \
       --add-host=database:$HOST_IP \
       --env-file foreman/environments/dev \
       dr_surveyor

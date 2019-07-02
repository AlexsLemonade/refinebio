#!/bin/bash

# Exit on error:
set -e

# Running this script will start an interactive python shell running
# within the context of a Docker container.
# By default the Docker container will be for the workers project.
# This can be changed by modifying the --env-file command line arg,
# changing workers/Dockerfile to the appropriate Dockerfile,
# and changing the volume_directory path.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd "$script_directory"

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/foreman/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
    chmod -R a+rwX "$volume_directory"
fi

docker build -t dr_shell -f foreman/dockerfiles/Dockerfile.foreman .

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

docker run -it \
       --add-host="database:$DB_HOST_IP" \
       --add-host="nomad:$HOST_IP" \
       --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
       --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
       --env-file foreman/environments/local \
       --volume /tmp:/tmp \
       --volume "$volume_directory":/home/user/data_store \
       --interactive dr_shell python3 manage.py shell

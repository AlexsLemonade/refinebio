#!/bin/sh

# Exit on error:
set -e

# Running this script will start an interactive python shell running
# within the context of a Docker container.
# By default the Docker container will be for the foreman project.
# This can be changed by modifying the --env-file command line arg,
# changing foreman/Dockerfile to the appropriate Dockerfile,
# and changing the volume_directory path.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Import functions in common.sh
. ./common.sh

# Get access to all of refinebio
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/../foreman/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

docker build \
    --file foreman/dockerfiles/Dockerfile.foreman \
    --tag dr_shell \
    .

DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
    --add-host="database:$DB_HOST_IP" \
    --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
    --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    --env-file foreman/environments/local \
    --interactive \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    --volume /tmp:/tmp \
    dr_shell \
    python3 manage.py shell

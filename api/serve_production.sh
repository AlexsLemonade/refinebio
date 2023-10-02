#!/bin/sh

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# However, in order to give Docker access to all the code we have to
# move up a level.
cd ..

./scripts/prepare_image.sh -i api -s api

. ./scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)
STATIC_VOLUMES=/tmp/volumes_static

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --detach \
    --env-file api/environments/local \
    --interactive \
    --link drdb:postgres \
    --platform linux/amd64 \
    --publish 8081:8081 \
    --tty \
    --volume "$STATIC_VOLUMES":/tmp/www/static \
    "$DOCKERHUB_REPO/dr_api" \
    /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"

#!/bin/sh

# Script for migrating the database using a Docker container so no
# virtual environment is needed on the host machine.

# Exit on error.
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

./prepare_image.sh -i base -s common
./prepare_image.sh -i migrations -s common

. ./common.sh

DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env-file ../common/environments/local \
    --interactive \
    --platform linux/amd64 \
    --volume "$script_directory/../common/data_refinery_common":/home/user/data_refinery_common \
    "$DOCKERHUB_REPO/dr_migrations:$SYSTEM_VERSION" \
    python3 manage.py makemigrations data_refinery_common

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env-file ../common/environments/local \
    --platform linux/amd64 \
    --volume "$script_directory/../common/data_refinery_common":/home/user/data_refinery_common \
    "$DOCKERHUB_REPO/dr_migrations:$SYSTEM_VERSION" \
    python3 manage.py migrate

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env-file ../common/environments/local \
    --platform linux/amd64 \
    --volume "$script_directory/../common/data_refinery_common":/home/user/data_refinery_common \
    "$DOCKERHUB_REPO/dr_migrations:$SYSTEM_VERSION" \
    python3 manage.py createcachetable

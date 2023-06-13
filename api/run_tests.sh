#!/bin/sh -e

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

# Ensure that Postgres is running.
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

# Ensure that ElasticSearch is running.
if ! [ "$(docker ps --filter name=dres -q)" ]; then
    echo "You must start elasticsearchfirst with:" >&2
    echo "./scripts/run_es.sh" >&2
    exit 1
fi

project_root=$(pwd) # "cd .." called above.
volume_directory="$project_root/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

./scripts/prepare_image.sh -i api_local -s api

. ./scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)
ES_HOST_IP=$(get_docker_es_ip_address)

# Only run interactively if we are on a TTY.
if [ -t 1 ]; then
    INTERACTIVE="--interactive"
fi

# shellcheck disable=SC2086
docker run \
    --add-host=database:"$DB_HOST_IP" \
    --add-host=elasticsearch:"$ES_HOST_IP" \
    --env-file api/environments/test \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    $INTERACTIVE \
    "$DOCKERHUB_REPO/dr_api_local" \
    bash -c "$(run_tests_with_coverage "$@")"

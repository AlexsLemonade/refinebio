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

./scripts/prepare_image.sh -i api_base -s api
./scripts/prepare_image.sh -i api_local -s api

. ./scripts/common.sh

# Only run interactively if we are on a TTY.
if [ -t 1 ]; then
    INTERACTIVE="--interactive"
fi

# shellcheck disable=SC2086
docker run \
    --network refinebio_default \
    --env-file api/environments/test \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    $INTERACTIVE \
    "$DOCKERHUB_REPO/dr_api_local:$SYSTEM_VERSION" \
    bash -c "$(run_tests_with_coverage "$@")"

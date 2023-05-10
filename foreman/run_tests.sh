#!/bin/sh

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# Set up the data volume directory if it does not already exist.
# Since the end-to-end tests are run from the Foreman image, use the
# top level test_volume rather than one nested within the foreman
# directory.
project_root=$(cd .. && pwd)
volume_directory="$project_root/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

# However, in order to give Docker access to all the code we have to
# move up a level.
cd ..

# First ensure Postgres is running.
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

./scripts/prepare_image.sh -i foreman -s foreman

. ./scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)

# Only run interactively if we are on a TTY.
if [ -t 1 ]; then
    INTERACTIVE="--interactive"
fi

# shellcheck disable=SC2086
docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env-file foreman/environments/test \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    $INTERACTIVE \
    "$DOCKERHUB_REPO/dr_foreman" \
    bash -c "$(run_tests_with_coverage --exclude-tag=manual "$@")"

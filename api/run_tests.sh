#!/bin/sh -e

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Ensure that postgres is running
if ! [ "$(docker ps --filter name=drdb -q)" ]; then
    echo "You must start Postgres first with:" >&2
    echo "./scripts/run_postgres.sh" >&2
    exit 1
fi

project_root=$(pwd) # "cd .." called above
volume_directory="$project_root/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
    chmod -R a+rwX "$volume_directory"
fi

./scripts/prepare_image.sh -i api_local -s api

. ./scripts/common.sh
HOST_IP=$(get_ip_address || echo 127.0.0.1)
DB_HOST_IP=$(get_docker_db_ip_address)
ES_HOST_IP=$(get_docker_es_ip_address)


# Only run interactively if we are on a TTY
if [ -t 1 ]; then
    INTERACTIVE="-i"
fi

docker run -t $INTERACTIVE \
       --add-host=database:"$DB_HOST_IP" \
       --add-host=nomad:"$HOST_IP" \
       --add-host=elasticsearch:"$ES_HOST_IP" \
       --env-file api/environments/test \
       --volume "$volume_directory":/home/user/data_store \
       ccdlstaging/dr_api_local bash -c "$(run_tests_with_coverage "$@")"

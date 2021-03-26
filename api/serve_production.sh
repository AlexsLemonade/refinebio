#!/bin/sh

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

./scripts/prepare_image.sh -i api_production -s api

. ./scripts/common.sh
DB_HOST_IP=$(get_docker_db_ip_address)

STATIC_VOLUMES=/tmp/volumes_static

docker run \
       --add-host=database:"$DB_HOST_IP" \
       --env-file api/environments/local \
       --link drdb:postgres \
       -v "$STATIC_VOLUMES":/tmp/www/static \
       -p 8081:8081 \
       -it -d ccdlstaging/dr_api_production /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"

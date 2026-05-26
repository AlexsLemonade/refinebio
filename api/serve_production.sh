#!/bin/sh

# Run the production API server locally using the uwsgi prod image
# (dr_api). Mirrors what runs on the production API EC2 box.
#
# Requires: rbio compose:up postgres
# (so the refinebio_default network and the `database` alias exist)

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

./bin/rbio build api

STATIC_VOLUMES=/tmp/volumes_static

docker run \
    --detach \
    --env-file api/environments/local \
    --interactive \
    --network refinebio_default \
    --platform linux/amd64 \
    --publish 8081:8081 \
    --tty \
    --volume "$STATIC_VOLUMES":/tmp/www/static \
    "$DOCKERHUB_REPO/dr_api" \
    /bin/sh -c "/home/user/collect_and_run_uwsgi.sh"

#!/bin/sh

# Script for running the Data Refinery Surveyor container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# However, in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

. ./scripts/common.sh

./scripts/prepare_image.sh -i base -s common
./scripts/prepare_image.sh -i foreman -s foreman

. ./scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
    --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    --env-file foreman/environments/local \
    --interactive \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    "$DOCKERHUB_REPO/dr_foreman:$SYSTEM_VERSION" \
    python3 manage.py "$@"

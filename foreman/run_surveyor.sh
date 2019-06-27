#!/bin/sh

# Script for running the Data Refinery Surveyor container

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
    chmod -R a+rwX "$volume_directory"
fi

./scripts/prepare_image.sh -i foreman -s foreman

. ./scripts/common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)
ES_HOST_IP=$(get_docker_es_ip_address)

docker run -it \
       --add-host=database:"$DB_HOST_IP" \
       --add-host=nomad:"$HOST_IP" \
       --add-host=elasticsearch:"$ES_HOST_IP" \
       --env-file foreman/environments/local \
       --env AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
       --env AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
       --volume "$volume_directory":/home/user/data_store \
       ccdlstaging/dr_foreman python3 manage.py "$@"

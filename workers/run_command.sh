#!/bin/sh

# Script for executing Django management commands within a Docker container.

# Exit on failure
set -e

while getopts "i:" opt; do
    case $opt in
        i)
            image=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

if [ -z "$image" ]; then
    image="smasher"
else
    shift
    shift
fi

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
    echo "./run_postgres.sh" >&2
    exit 1
fi

volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
    chmod -R a+rwX "$volume_directory"
fi

. ./common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

chmod -R a+rwX "$volume_directory"

./prepare_image.sh -i "$image" -s workers
image_name=ccdlstaging/dr_"$image"

docker run \
       --add-host=database:"$DB_HOST_IP" \
       --add-host=nomad:"$HOST_IP" \
       --env-file workers/environments/local \
       --env AWS_ACCESS_KEY_ID \
       --env AWS_SECRET_ACCESS_KEY \
       --volume "$volume_directory":/home/user/data_store \
       --link drdb:postgres \
       -it "$image_name" Rscript -e "source('more_accessions.R'); get_random_sample_accessions('SRA_supported_runs_human_mouse_rat/mouse.xml', '/home/user/data_store/more_mouse.txt', blacklist='crunch_lists/mouse_random_0.01_SRA_study_accessions.txt')"

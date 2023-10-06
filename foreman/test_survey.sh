#!/bin/sh

# Script for testing the surveying of an accession manually (e.g. from a dataset request)
print_options() {
    cat <<EOF
There are two required arguments for this script that come from survey_experiment():
-s specifies which surveyor to use.
-a specifies which experiment accession should be surveyed.
EOF
}

while getopts ":s:a:h" opt; do
    case $opt in
    s)
        export SURVEYOR="$OPTARG"
        ;;
    a)
        export ACCESSION="$OPTARG"
        ;;
    h)
        print_description
        echo
        print_options
        exit 0
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        print_options >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        print_options >&2
        exit 1
        ;;
    esac
done

if [ -z "$SURVEYOR" ]; then
    echo 'Error: must specify which surveyor to use with -s' >&2
    exit 1
fi

if [ -z "$ACCESSION" ]; then
    echo 'Error: must specify which accession to survey with -a' >&2
    exit 1
fi

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

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --env ACCESSION="$ACCESSION" \
    --env SURVEYOR="$SURVEYOR" \
    --env-file foreman/environments/test \
    --interactive \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    "$DOCKERHUB_REPO/dr_foreman" \
    bash -c "python3 manage.py test --tag=manual ."

#!/bin/sh

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Set up the data volume directory if it does not already exist.
# Since the end-to-end tests are run from the Foreman image, use the
# top level test_volume rather than one nested within the foreman
# directory.
project_root=$(cd .. && pwd)
volume_directory="$project_root/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
    chmod -R a+rwX "$volume_directory"
fi

# temp for testing locally.
../scripts/prepare_image.sh -i foreman -s foreman


while read -r row; do
    # Exporting an expansion rather than a variable, which is exactly what we want to do.
    # shellcheck disable=SC2163
    export "${row}"
done < ../infrastructure/prod_env

# Hardcode ccdlstaging because this should only ever be run in staging.
docker run -t \
       --env-file ../infrastructure/prod_env \
       --env RUNNING_IN_CLOUD=False \
       --env DATABASE_HOST="$DATABASE_PUBLIC_HOST" \
       --env JOB_DEFINITION_PREFIX="$USER_$STAGE_" \
       --env REFINEBIO_BASE_URL="http://$API_HOST/v1/" \
       --volume "$volume_directory":/home/user/data_store \
       ccdlstaging/dr_foreman python3 manage.py test --no-input --parallel=2 --testrunner='data_refinery_foreman.test_runner.NoDbTestRunner' data_refinery_foreman.foreman.test_end_to_end
